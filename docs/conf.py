# type: ignore
# pylint: skip-file

# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
from os import path as op
import sys

sys.path.insert(0, op.abspath(op.join("..")))
print(sys.path)

# -- Project information -----------------------------------------------------

project = "symplyphysics"
copyright = "2024, Symplyphysics"
author = "blackyblack"

# The full version, including alpha/beta/rc tags
release = "1.0.0"

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named "sphinx.ext.*") or your custom
# ones.
extensions = [
    "sphinx.ext.viewcode",
    "sphinx.ext.autosectionlabel",
    "sphinx.ext.mathjax",
    "sphinx.ext.githubpages",
    "sphinx_sitemap",
]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "alabaster"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]

html_extra_path = ["_extra/googlea3669135f4e8d2d0.html"]

# If true, the current module name will be prepended to all description
# unit titles (such as .. function::).
add_module_names = False

mathjax_path = "scipy-mathjax/MathJax.js?config=scipy-mathjax"

# -- Options for sphinx-sitemap

html_baseurl = "https://symplyphysics.github.io/"
sitemap_url_scheme = "{link}"
