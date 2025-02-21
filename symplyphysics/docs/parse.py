import ast
import re
from dataclasses import dataclass
from enum import Enum
from typing import Any, Optional
from sympy.physics.units.systems.si import SI
from ..core.symbols.symbols import DimensionSymbol, Function
from .printer_code import code_str
from .printer_latex import latex_str


class LawDirectiveTypes(Enum):
    SYMBOL = 0
    LATEX = 1


class LawSymbolTypes(Enum):
    SYMBOL = 0
    FUNCTION = 1


@dataclass
class LawDirective:
    start: int
    end: int
    directive_type: LawDirectiveTypes


@dataclass
class LawSymbol:
    symbol: str
    symbol_type: LawSymbolTypes
    latex: Optional[str]
    dimension: str


@dataclass
class MemberWithDoc:
    name: str
    docstring: str
    symbol: Optional[LawSymbol]
    directives: list[LawDirective]
    value: Any


@dataclass
class FunctionWithDoc:
    name: str
    parameters: list[str]
    returns: Optional[str]
    docstring: str


_LAWS_SYMPY_EVAL_PATTERN = re.compile(r"\n?:laws:sympy-eval::\n?")


def _docstring_clean(doc: str) -> str:
    doc = _LAWS_SYMPY_EVAL_PATTERN.sub("", doc)
    doc = doc.strip("\n")
    return doc


_LAWS_SYMBOL_STR = ":laws:symbol::"
_LAWS_LATEX_STR = ":laws:latex::"


def _docstring_find_law_directives(doc: str) -> list[LawDirective]:
    directives: list[LawDirective] = []

    position = doc.find(_LAWS_SYMBOL_STR)
    if position >= 0:
        directives.append(
            LawDirective(position, position + len(_LAWS_SYMBOL_STR), LawDirectiveTypes.SYMBOL))

    position = doc.find(_LAWS_LATEX_STR)
    if position >= 0:
        directives.append(
            LawDirective(position, position + len(_LAWS_LATEX_STR), LawDirectiveTypes.LATEX))

    return directives


def find_title_and_description(doc: str) -> Optional[tuple[str, str]]:

    def is_section_break(line: str, char: str) -> bool:
        return all(c == char for c in line)

    lines = doc.splitlines()

    for idx, line in enumerate(lines[1:], start=1):
        if not line:
            continue

        if is_section_break(line, "=") or is_section_break(line, "-"):
            section_break = idx
            break
    else:
        return None

    title = "\n".join(lines[:section_break])
    description_lines = lines[section_break + 1:]

    # Remove possible empty lines after section break
    while description_lines:
        if description_lines[0]:
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
        doc = _docstring_clean(doc)

        name = stmt.name
        parameters = [arg.arg for arg in stmt.args.args]

        returns = stmt.returns
        return_type = returns.id if returns and isinstance(returns, ast.Name) else None

        return FunctionWithDoc(name, parameters, return_type, doc)

    def process_assign(stmt: ast.Assign) -> Optional[str]:
        for target in stmt.targets:
            if isinstance(target, ast.Name):
                name = target.id
                break
        else:
            return None

        return name

    def process_body() -> tuple[list[FunctionWithDoc], list[str], dict[str, str]]:
        functions: list[FunctionWithDoc] = []
        member_names: list[str] = []
        current_member_name: Optional[str] = None
        docstrings: dict[str, str] = {}

        for e in module.body:
            if isinstance(e, ast.FunctionDef):
                function_ = process_function(e)
                if function_:
                    functions.append(function_)
                continue

            if isinstance(e, ast.Assign):
                current_member_name = process_assign(e)
                if current_member_name:
                    member_names.append(current_member_name)
                continue

            if (isinstance(e, ast.Expr) and current_member_name and
                    isinstance(e.value, ast.Constant)):
                docstrings[current_member_name] = e.value.value

        return functions, member_names, docstrings

    def process_member_name(name: str) -> Optional[MemberWithDoc]:
        doc = docstrings.get(name)
        if doc is None:
            return None
        doc = _docstring_clean(doc)

        value = context[name]
        symbol: Optional[LawSymbol] = None

        if isinstance(value, DimensionSymbol):
            # TODO: use `core.dimensions.print_dimension` here when merged
            dimension = "dimensionless" if SI.get_dimension_system().is_dimensionless(
                value.dimension) else str(value.dimension.name)

            symbol_name = code_str(value)
            symbol_latex = latex_str(value)

            symbol_type = (LawSymbolTypes.FUNCTION
                if isinstance(value, Function) else LawSymbolTypes.SYMBOL)

            symbol = LawSymbol(symbol_name, symbol_type, symbol_latex, dimension)

        directives = _docstring_find_law_directives(doc)
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
