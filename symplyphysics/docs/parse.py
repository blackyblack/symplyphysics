"""
This module provides a parsing functionality for the documentation of laws.

Classes:

* `LawDirective` contains the location information of a directive in the docstring.
  `LawDirectiveType` enumerates possible types of directives.

* `LawSymbol` stores the printing information and dimension of a module's variable.
  `LawSymbolType` enumerates possible types of symbols.

* `MemberWithDoc` contains the documentation information as well as the associated `LawSymbol`.

* `FunctionWithDoc` contains the documentation information and the function signature.

Functions:

* `find_title_and_description` attempts to find a title and a description in the top docstring of
  a module.

* `find_members_and_functions` accepts an `ast.Module` object and locates the documented variables
  and functions within it.
"""

import ast
import re
from dataclasses import dataclass
from enum import Enum
from typing import Any, Optional
from sympy.physics.units.systems.si import SI
from ..core.symbols.symbols import DimensionSymbol, Function
from .printer_code import code_str
from .printer_latex import latex_str


class LawDirectiveType(Enum):
    """Enumeration of possible directives."""

    SYMBOL = 0
    """Code-printed representation."""

    LATEX = 1
    """Latex-printed representation."""


class LawSymbolType(Enum):
    """Enumeration of possible types of module variables."""

    SYMBOL = 0
    """Corresponds to the `symplyphysics.Symbol` and `symplyphysics.IndexedSymbol` classes."""

    FUNCTION = 1
    """Corresponds to the `symplyphysics.Function` class."""


@dataclass
class LawDirective:
    """Represents the location information within the docstring for a directive."""

    start: int
    """Initial location of the directive."""

    end: int
    """Final location of the directive."""

    directive_type: LawDirectiveType
    """Type of the directive."""


@dataclass
class LawSymbol:
    """Represents the symbol or function associated with a module variable."""

    symbol: str
    """Code-printed representation of the symbol."""

    symbol_type: LawSymbolType
    """Type of the symbol."""

    latex: Optional[str]
    """Latex-printed representation of the symbol."""

    dimension: str
    """Dimension of the symbol."""


@dataclass
class MemberWithDoc:
    """Represents a module variable."""

    name: str
    """Name of the variable."""

    docstring: str
    """Docstring attached to the variable."""

    symbol: Optional[LawSymbol]
    """The symbol information associated with the variable."""

    directives: list[LawDirective]
    """Directived parsed from the documentation."""

    value: Any
    """The actual symbol associated with the variable obtained via `exec`."""


@dataclass
class FunctionWithDoc:
    """Represents a function defined in the module."""

    name: str
    """Name of the function."""

    parameters: list[str]
    """Parameter list of the function."""

    returns: Optional[str]
    """Return type of the function if it is a simple type, else `None`."""

    docstring: str
    """Docstring attached to the function."""


_LAWS_SYMPY_EVAL_PATTERN = re.compile(r"\n?:laws:sympy-eval::\n?")


def _clean_docstring(doc: str) -> str:
    doc = _LAWS_SYMPY_EVAL_PATTERN.sub("", doc)
    doc = doc.strip("\n")
    return doc


_LAWS_SYMBOL_STR = ":laws:symbol::"
_LAWS_LATEX_STR = ":laws:latex::"


def _find_law_directives(doc: str) -> list[LawDirective]:
    """
    Attempts to locate the law directives within ``doc``.
    
    Returns a list of 0, 1, or 2 directives.
    """

    directives: list[LawDirective] = []

    position = doc.find(_LAWS_SYMBOL_STR)
    if position >= 0:
        directives.append(
            LawDirective(position, position + len(_LAWS_SYMBOL_STR), LawDirectiveType.SYMBOL))

    position = doc.find(_LAWS_LATEX_STR)
    if position >= 0:
        directives.append(
            LawDirective(position, position + len(_LAWS_LATEX_STR), LawDirectiveType.LATEX))

    return directives


def find_title_and_description(doc: str) -> Optional[tuple[str, str]]:
    """
    Locates the title and description of ``doc``.
    
    Returns `None` if no section break is present.
    """

    def is_section_break(line: str, char: str) -> bool:
        return all(c == char for c in line)

    lines = doc.splitlines()

    section_break: Optional[int] = None

    for idx, line in enumerate(lines[1:], start=1):
        if not line:
            continue

        if is_section_break(line, "=") or is_section_break(line, "-"):
            section_break = idx
            break

    if section_break is None:
        return None

    title = "\n".join(lines[:section_break])
    description_lines = lines[section_break + 1:]

    # Remove possible empty lines after section break
    while True:
        if not description_lines or description_lines[0]:
            break
        description_lines.pop(0)

    description = "\n".join(description_lines)

    return title, description


def find_members_and_functions(
    module: ast.Module,) -> tuple[list[MemberWithDoc], list[FunctionWithDoc]]:
    """
    Parses ``module`` collecting documented members (i.e. variables) and defined functions.
    """

    def process_function(stmt: ast.FunctionDef) -> Optional[FunctionWithDoc]:
        doc = ast.get_docstring(stmt)
        if doc is None:
            return None
        doc = _clean_docstring(doc)

        name = stmt.name
        parameters = [arg.arg for arg in stmt.args.args]

        returns = stmt.returns
        return_type = returns.id if returns and isinstance(returns, ast.Name) else None

        return FunctionWithDoc(name, parameters, return_type, doc)

    def process_assign(stmt: ast.Assign) -> Optional[str]:
        for target in stmt.targets:
            if isinstance(target, ast.Name):
                return target.id

        return None

    def process_body() -> tuple[list[FunctionWithDoc], list[str], dict[str, str]]:
        functions: list[FunctionWithDoc] = []
        member_names: list[str] = []
        current_member_name: Optional[str] = None
        docstrings: dict[str, str] = {}

        for stmt in module.body:
            if isinstance(stmt, ast.FunctionDef):
                function_ = process_function(stmt)
                if function_:
                    functions.append(function_)
                continue

            if isinstance(stmt, ast.Assign):
                current_member_name = process_assign(stmt)
                if current_member_name:
                    member_names.append(current_member_name)
                continue

            if (isinstance(stmt, ast.Expr) and current_member_name and
                    isinstance(stmt.value, ast.Constant)):
                docstrings[current_member_name] = stmt.value.value

        return functions, member_names, docstrings

    def process_member_name(name: str) -> Optional[MemberWithDoc]:
        doc = docstrings.get(name)
        if doc is None:
            return None
        doc = _clean_docstring(doc)

        value = context[name]
        symbol: Optional[LawSymbol] = None

        if isinstance(value, DimensionSymbol):
            # TODO: use `core.dimensions.print_dimension` here when merged
            dimension = "dimensionless" if SI.get_dimension_system().is_dimensionless(
                value.dimension) else str(value.dimension.name)

            symbol_name = code_str(value)
            symbol_latex = latex_str(value)

            symbol_type = (LawSymbolType.FUNCTION
                if isinstance(value, Function) else LawSymbolType.SYMBOL)

            symbol = LawSymbol(symbol_name, symbol_type, symbol_latex, dimension)

        directives = _find_law_directives(doc)
        return MemberWithDoc(name, doc, symbol, directives, value)

    functions, member_names, docstrings = process_body()

    compiled = compile(module, "string", "exec")
    context: dict[str, Any] = {}
    exec(compiled, {}, context)  # pylint: disable=exec-used

    members: list[MemberWithDoc] = []
    for member_name in member_names:
        member = process_member_name(member_name)
        if member:
            members.append(member)

    return members, functions


__all__ = [
    "LawDirectiveType",
    "LawSymbolType",
    "LawDirective",
    "LawSymbol",
    "MemberWithDoc",
    "FunctionWithDoc",
    "find_title_and_description",
    "find_members_and_functions",
]
