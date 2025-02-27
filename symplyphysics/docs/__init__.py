"""
This package provides the functionality for printing code and building documentation.
"""

from .printer_code import code_str
from .printer_latex import latex_str

__all__ = [
    # only code-printing functions are exported
    "code_str",
    "latex_str",
]
