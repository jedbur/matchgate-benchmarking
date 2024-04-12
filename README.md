# matchgate-benchmarking

This is supplementary code for the paper [*"A Lightweight Protocol for Matchgate Fidelity Estimation"*](https://arxiv.org/abs/2404.07974), used to simulate Algorithm 1 as described in the text, using Python 3.11 and qiskit. We recommend using a virtual environment to install the project, however the code can be directly copied from the `mgbenchmark` folder.

## Project Structure

The source code is located in the `mgbenchmark` folder, which contains the `__init__.py`, `main.py`, and `utils.py` files.

The `examples` folder contains `.csv` files with quoted data, and four Jupyter Notebooks:

* `matchgate_benchmarking.ipynb`: A demo of the Matchgate Fidelity Estimation protocol (with brief discussion).
* `matchgate_benchmarking_data_collection.ipynb`: A demo of the protocol with different noise models, as well as code used to collect data for the paper.
* `matchgate_tomography.ipynb`: A demo of efficient matchgate tomography, first described by [Oszmaniec, Dangniam, Morales and Zimboras](https://arxiv.org/abs/2012.15825).
* `data_plots.ipynb`: Code used to make the figures.

The project uses toml for configuration instead of setup.py. The configuration file is located in `pyproject.toml`. We've used Pyright for static type checking, pre-commit for code formatting, Black for code formatting and Ruff for linting. The configuration for these tools is located in the `.pre-commit-config.yaml` and `ruff.yaml` file. This was included for convenience, and can be ignored.

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
