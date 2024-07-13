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
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.autosectionlabel",
    "sphinx.ext.autosummary",
    "sphinx.ext.mathjax",
    "sphinx.ext.githubpages",
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ["templates"]

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

# If true, the current module name will be prepended to all description
# unit titles (such as .. function::).
add_module_names = False

napoleon_custom_sections = ["Symbol", "Symbols", "Latex"]

autosummary_generate = True

# Monkey-patch autosummary template context
from sphinx.ext.autosummary.generate import AutosummaryRenderer


def get_first_line(docstring: str | None) -> str:
    if not docstring:
        # __doc__ can be None
        return ""
    lines = docstring.split("\n")
    return lines[0]


def extract_module_docstring(mod_name) -> str:
    """See _templates/autosummary/base.rst"""
    import sys
    mod = sys.modules[mod_name]
    return get_first_line(getattr(mod, "__doc__", ""))


def partial_name(fullname):
    parts = fullname.split(".")
    return parts[-1]


# Patch autosummary internals to allow our tuned templates to access
# necessary Python functions
def fixed_init(self, app):
    AutosummaryRenderer.__old_init__(self, app)
    self.env.filters["extract_module_docstring"] = extract_module_docstring
    self.env.filters["partial_name"] = partial_name


AutosummaryRenderer.__old_init__ = AutosummaryRenderer.__init__
AutosummaryRenderer.__init__ = fixed_init

mathjax_path = "scipy-mathjax/MathJax.js?config=scipy-mathjax"
