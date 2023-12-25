# matchgate-benchmarking

This is code for matchgate fidelity estimation, using Python 3.11 and qiskit. It implements the algorithm described in our upcoming paper. We recommend using a virtual environment to install the project, however the code can be directly copied from the `mgbenchmark` folder.

## Project Structure

The source code is located in the `mgbenchmark` folder, which contains the `__init__.py`, `main.py`, and `utils.py` files.

The algorithm and discussion is located in the `matchgate_benchmarking` Jupyter Notebook, which is located in the `examples` folder.

The project uses toml for configuration instead of setup.py. The configuration file is located in `pyproject.toml`.

The project includes Pyright for static type checking, pre-commit for code formatting, Black for code formatting and Ruff for linting. The configuration for these tools is located in the `.pre-commit-config.yaml` and `ruff.yaml` file. This is included for my own convenience when writing, and can be ignored.

## Installation

To install the project, clone the repository and run:

```sh
python -m venv .venv
source .venv/bin/activate
pip install -U pip setuptools
pip install -r requirements.txt
pre-commit install
```

Then install the project using:

```sh
pip install -e .
```
