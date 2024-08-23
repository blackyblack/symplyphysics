import ast
from dataclasses import dataclass
from enum import Enum
import importlib
from typing import Any, Optional
from sympy.physics.units.systems.si import SI

from ..core.symbols.symbols import DimensionSymbol


class LawDirectiveTypes(Enum):
    SYMBOL = 0
    LATEX = 1

@dataclass
class LawDirective:
    start: int
    end: int
    directive_type: LawDirectiveTypes

@dataclass
class LawSymbol:
    symbol: str
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


def _docstring_clean(doc: str) -> str:
    doc = doc.replace("\r\n", "\n")
    doc = doc.replace("\n\r", "\n")
    while True:
        if not doc.startswith("\n"):
            break
        doc = doc.removeprefix("\n")
    while True:
        if not doc.endswith("\n"):
            break
        doc = doc.removesuffix("\n")
    return doc


def _docstring_find_law_directives(doc: str) -> list[LawDirective]:
    directives = []
    position = doc.find(":laws:symbol::\n")
    if position >= 0:
        directives.append(LawDirective(position, position + len(":laws:symbol::\n"), LawDirectiveTypes.SYMBOL))
    position = doc.find(":laws:latex::\n")
    if position >= 0:
        directives.append(LawDirective(position, position + len(":laws:latex::\n"), LawDirectiveTypes.LATEX))
    return directives


def find_title(content: str) -> Optional[str]:
    content_lines = content.splitlines()
    for i, c in enumerate(content_lines):
        l = len(c)
        if l == 0:
            continue
        if c == ("=" * l) and i > 0:
            return "\n".join(content_lines[0:i])
        if c == ("-" * l) and i > 0:
            return "\n".join(content_lines[0:i])
    return None


def find_description(content: str) -> Optional[str]:
    content_lines = content.splitlines()
    section_break = 0
    for i, c in enumerate(content_lines):
        l = len(c)
        if l == 0:
            continue
        if c == ("=" * l) and i > 0:
            section_break = i
            break
        if c == ("-" * l) and i > 0:
            section_break = i
            break
    if section_break == 0:
        return None
    content_lines = content_lines[section_break + 1:]
    # Remove possible empty lines after title separator
    while True:
        if len(content_lines) == 0:
            return ""
        if len(content_lines[0]) == 0:
            content_lines.pop(0)
            continue
        break
    return "\n".join(content_lines)


def find_members_and_functions(module_name: str) -> list[MemberWithDoc | FunctionWithDoc]:
    law_functions: list[FunctionWithDoc] = []
    law_members: list[str] = []
    current_member: Optional[str] = None
    docstrings: dict[str, str] = {}
    module = importlib.import_module(module_name)
    if module.__file__ is None:
        print(f"Unable to read {module} contents")
        return []
    with open(module.__file__, "r", encoding="utf-8") as f:
        package_content = f.read()
    content = ast.parse(package_content)
    for e in content.body:
        if isinstance(e, ast.FunctionDef):
            doc = ast.get_docstring(e)
            if doc is None:
                continue
            args_list_str: list[str] = [a.arg for a in e.args.args]
            maybe_return_str = None if e.returns is None else (e.returns.id if isinstance(e.returns, ast.Name) else None)
            doc = _docstring_clean(doc)
            law_functions.append(FunctionWithDoc(e.name, args_list_str, maybe_return_str, doc))
            continue
        if isinstance(e, ast.Assign):
            for t in e.targets:
                try:
                    name = getattr(t, "id")
                except AttributeError:
                    continue
                law_members.append(name)
                current_member = name
            continue
        if isinstance(e, ast.Expr) and current_member is not None and isinstance(
                e.value, ast.Constant):
            docstrings[current_member] = e.value.value
            continue
    result: list[MemberWithDoc | FunctionWithDoc] = []
    for v in law_members:
        doc = docstrings.get(v)
        if doc is None:
            continue
        doc = _docstring_clean(doc)
        sym = getattr(module, v)
        law_symbol = None
        if isinstance(sym, DimensionSymbol):
            dimension = "dimensionless" if SI.get_dimension_system().is_dimensionless(sym.dimension) else str(sym.dimension.name)
            law_symbol = LawSymbol(sym.display_symbol, sym.display_latex, dimension)
        directives = _docstring_find_law_directives(doc)
        result.append(MemberWithDoc(v, doc, law_symbol, directives, sym))
    for lf in law_functions:
        result.append(lf)
    return result
