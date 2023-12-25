"""Utility functions for mgbenchmark."""

from itertools import combinations

import numpy as np
from numpy.typing import NDArray
from qiskit.quantum_info import Pauli  # type: ignore[import]


def generate_jw_list(n: int) -> list[NDArray[np.complex128]]:
    """Returns the Jordan-Wigner generating set for n qubits.

    Args:
        n: The number of qubits.

    Returns:
        jw_list: A list of 2n matrices, dimension 2n x 2n.
    """
    pauli_id = np.array([[1, 0], [0, 1]], dtype=np.complex128)
    pauli_x = np.array([[0, 1], [1, 0]], dtype=np.complex128)
    pauli_y = np.array([[0, -1j], [1j, 0]], dtype=np.complex128)
    pauli_z = np.array([[1, 0], [0, -1]], dtype=np.complex128)
    jw_list: list[NDArray[np.complex128]] = []

    for k in range(1, n + 1, 1):
        c_1 = 1
        c_2 = 1
        for _j in range(k - 1):
            c_1 = np.kron(c_1, pauli_z)
            c_2 = c_1.copy()
        c_1 = np.kron(c_1, pauli_x)
        c_2 = np.kron(c_2, pauli_y)
        for _j in range(n - k):
            c_1 = np.kron(c_1, pauli_id)
            c_2 = np.kron(c_2, pauli_id)
        jw_list.append(c_1)
        jw_list.append(c_2)

    return jw_list


def generate_jw_list_pauli(n: int) -> list[Pauli]:
    """Returns the Jordan-Wigner generating set for n qubits.

    Args:
        n: The number of qubits.

    Returns:
        jw_list: A list of 2n Pauli Operators, dimension 2n x 2n.
    """
    jw_pauli_list: list[Pauli] = []

    for k in range(1, n + 1, 1):
        c_1 = ""
        c_2 = ""
        for _j in range(k - 1):
            c_1 = c_1 + "Z"
            c_2 = c_2 + "Z"
        c_1 = c_1 + "X"
        c_2 = c_2 + "Y"
        for _j in range(n - k):
            c_1 = c_1 + "I"
            c_2 = c_2 + "I"
        jw_pauli_list.append(Pauli(c_1))
        jw_pauli_list.append(Pauli(c_2))

    return jw_pauli_list


def generate_jw_basis(n: int) -> list[NDArray[np.complex128]]:
    """Returns the full Jordan-Wigner basis for n qubits.

    Args:
        n: The number of qubits.

    Returns:
        jw_list: A list of 2^2n matrices, dimension 2n x 2n.
    """
    # Initialize the Jordan-Wigner generators
    jw_list = generate_jw_list(n)

    # Generate the remaining basis vectors
    c_0 = np.eye(2**n, dtype=np.complex128)
    jw_list = [c_0, *jw_list]

    for i in range(2, 2 * n + 2, 1):
        for s in generate_binary_strings(2 * n, i):
            c_s = np.eye(2**n, dtype=np.complex128)
            for b in range(2 * n):
                if s[b] == "1":
                    c_s = c_s @ jw_list[b + 1]
                else:
                    pass
            jw_list.append(c_s)

    return jw_list


def generate_binary_strings(n: int, k: int) -> list[str]:
    """Generates all binary strings of length n with Hamming weight k.

    Args:
        n: The length of the binary strings.
        k: The Hamming weight of the binary strings.

    Returns:
        nk_list: A list of (n, k)-binary strings, in lexographical order.
    """
    nk_list: list[str] = []
    for indices in combinations(range(n), k):
        # Create a binary string with '0's
        binary_string = ["0"] * n

        # Set the positions in the combination to '1'
        for index in indices:
            binary_string[index] = "1"

        # Convert the list of characters back to a string
        nk_list.append("".join(binary_string))

    return nk_list


def compound_matrix(u: NDArray[np.complex128], r: int) -> NDArray[np.complex128]:
    """Returns the k^th compound matrix of u.

    Compound matrices are defined as https://en.wikipedia.org/wiki/Compound_matrix.
    The columns I, J are subsets of {1, ..., n} of size k. The entries c_{I, J} are
    the determinants of the submatrices of u formed by preserving the rows of u indexed
    by I and the columns of u indexed by J.

    Args:
        u: An n x n square matrix.
        r: The size of the subsets.

    Returns:
        c: An (n choose k) x (n choose k) square matrix.
    """
    if u.shape[0] != u.shape[1]:
        msg = "R must be square."
        raise ValueError(msg)
    if r < 0 or r > u.shape[0]:
        return np.zeros((1, 1), dtype=np.complex128)

    n = u.shape[0]

    # Generate the list of all subsets of {1, ..., n} of size k, using binary strings
    subsets = generate_binary_strings(n, r)
    dim = len(subsets)

    # Convert the binary strings to lists of indices
    c = np.zeros((dim, dim), dtype=np.complex128)
    for i in subsets:
        for j in subsets:
            x: list[int] = []
            y: list[int] = []
            for k in range(0, n, 1):
                if i[k] == "0":
                    x = [*x, k]
                if j[k] == "0":
                    y = [*y, k]
            # Calculate the submatrix determinants
            mat = u.copy()
            mat = np.delete(mat, x, axis=0)
            mat = np.delete(mat, y, axis=1)
            c[subsets.index(i), subsets.index(j)] = np.linalg.det(mat)

    return c


def superop_to_dictionaries(
    superop: NDArray[np.complex128],
) -> tuple[dict[tuple[str, str], float], dict[tuple[str, str], float]]:
    """Converts a superoperator to dictionaries used in the benchmark protocol.

    Args:
        superop: The superoperator to be converted.

    Returns:
        dictionary_abs: A dictionary of the absolute values of non-zero coefficients of
        the superoperator, indexed by tuples of binary strings.
        dictionary_ph: A dictionary of the phases of non-zero coefficients.
        dictionary_prob: A dictionary of the probabilities of non-zero coefficients,
        calculated from the matrix elements.
    """
    n = int(np.log2(superop.shape[0]) // 2)

    bits = []
    for i in range(2 * n + 1):
        bits: list[str] = bits + generate_binary_strings(2 * n, i)

    dictionary_mat: dict[tuple[str, str], float] = {}
    dictionary_prob: dict[tuple[str, str], float] = {}
    for i in bits:
        for j in bits:
            if round(abs(superop[bits.index(i), bits.index(j)]), 15) != 0:
                dictionary_mat[(i, j)] = superop[bits.index(i), bits.index(j)]
                dictionary_prob[(i, j)] = (
                    abs(superop[bits.index(i), bits.index(j)]) ** 2
                ) / (2 ** (2 * n))

    return dictionary_mat, dictionary_prob


def dictionary_to_distribution(
    dictionary_prob: dict[tuple[str, str], float],
) -> tuple[list[tuple[str, str]], list[float]]:
    """Converts a dictionary of probabilities to a distribution.

    Args:
        dictionary_prob: A dictionary of probabilities, indexed by tuples of bitstrings.

    Returns:
        keys: A list of the keys of the dictionary.
        probabilities: A list of the probabilities of the keys.
    """
    keys = list(dictionary_prob.keys())
    probabilities = list(dictionary_prob.values())

    if round(sum(probabilities), 5) != 1:
        msg = "Probabilities must sum to 1"
        raise ValueError(msg)

    return keys, probabilities


def string_to_pauli(jw_list_pauli: list[Pauli], bitstring: str) -> Pauli:
    """Converts a binary string to a Pauli operator and its phase (multiple of pi/2)."""
    n = len(bitstring) // 2
    if len(bitstring) != len(jw_list_pauli):
        msg = "Mismatched bitstring and jw_list lengths."
        raise ValueError(msg)

    s_pauli = Pauli("I" * n)

    for i in range(2 * n):
        if bitstring[i] == "1":
            s_pauli: Pauli = jw_list_pauli[i].compose(s_pauli)  # type: ignore[assignment]

    return s_pauli
