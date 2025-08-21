# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

from pathlib import Path
import sys, os

ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(ROOT))

project = "Computable Model of Chemistry Textbook"
copyright = "2025, Peter G. Chang"
author = "Peter G. Chang"

version = "0.1.0"
release = "0.1.0"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "myst_parser",  # Markdown support
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.intersphinx",
    "sphinx.ext.viewcode",
    "sphinx_copybutton",
    "sphinx_autodoc_typehints",
    "sphinx.ext.napoleon",  # Google/NumPy docstrings
]

html_theme = "furo"

# Autodoc / typing / summaries
autosummary_generate = True
autodoc_typehints = "description"
autodoc_default_options = {
    "members": True,
    "undoc-members": True,
    "show-inheritance": True,
}

# If heavy deps cause import errors during autodoc, mock them:
autodoc_mock_imports = ["rdkit", "torch"]

# Docstring style
napoleon_google_docstring = True
napoleon_numpy_docstring = True

# Cross-links
intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    # "numpy": ("https://numpy.org/doc/stable/", None),
}