"""This is the main module of mgbenchmark."""
import numpy as np
import scipy  # type: ignore[import]
from numpy.typing import NDArray

from mgbenchmark.utils import (
    compound_matrix,
    generate_jw_basis,
    generate_jw_list,
)


def unitary_to_superoperator(u: NDArray[np.complex128]) -> NDArray[np.complex128]:
    """Converts a unitary to its superoperator matrix.

    Uses the Jordan-Wigner representation of the clifford algebra.

    Args:
        u: A 2^n x 2^n Unitary matrix.

    Returns:
        super_u: A 2^(2n) x 2^(2n) superoperator matrix.
    """
    n = int(np.log2(len(u)))

    # Initialize the Jordan-Wigner Basis
    jw_list = generate_jw_basis(n)

    # Calculate the matrix elements
    r = np.zeros((2 ** (2 * n), 2 ** (2 * n)), dtype=np.complex128)
    for i in range(2 ** (2 * n)):
        for j in range(2 ** (2 * n)):
            c_i = jw_list[i]
            c_j = jw_list[j]
            r[i, j] = np.trace(c_i.conj().T @ u @ c_j @ u.conj().T) / (2**n)

    return r


def mg_unitary_to_so(u: NDArray[np.complex128]) -> NDArray[np.complex128]:
    """Converts a matchgate unitary U to its Special Orthogonal matrix R.

    Uses the Jordan-Wigner representation of the clifford algebra.

    Args:
        u: A 2^n x 2^n Unitary matrix.

    Returns:
        r: A 2n x 2n Special Orthogonal matrix.
    """
    n = int(np.log2(len(u)))

    # Initialize the Jordan-Wigner Basis
    jw_list = generate_jw_list(n)

    # Calculate the SO matrix elements
    r = np.zeros((2 * n, 2 * n), dtype=np.complex128)
    for i in range(2 * n):
        for j in range(2 * n):
            c_i = jw_list[i]
            c_j = jw_list[j]
            r[i, j] = np.trace(c_i.conj().T @ u @ c_j.conj().T @ u.conj().T) / (2**n)

    return r


def mg_so_to_superoperator(r: NDArray[np.complex128]) -> NDArray[np.complex128]:
    """Converts a matchgate SO(2n) matrix R to its superoperator matrix.

    Uses compound matrix methods; should be more efficient to unitary_to_superoperator.

    Args:
        r: A 2n x 2n Special Orthogonal matrix.

    Returns:
        super_r: A 2^(2n) x 2^(2n) superoperator matrix.
    """
    if r.shape[0] != r.shape[1]:
        msg = "R must be square."
        raise ValueError(msg)

    if r.shape[0] % 2 != 0:
        msg = "R must be even dimensional."
        raise ValueError(msg)

    if (r @ r.T != np.eye(r.shape[0])).all():
        msg = "R must be orthogonal."
        raise ValueError(msg)

    n = r.shape[0] // 2

    super_r = scipy.linalg.block_diag([1], r)

    for k in range(2, 2 * n, 1):
        c = compound_matrix(r, k)
        super_r = scipy.linalg.block_diag(super_r, c)

    return scipy.linalg.block_diag(super_r, [1])
