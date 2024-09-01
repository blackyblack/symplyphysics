from typing import Sequence
from .printer_latex import latex_str
from .printer_code import code_str
from .parse import FunctionWithDoc, LawDirectiveTypes, MemberWithDoc


INDENT_SPACES = 4


def _indent_docstring(doc: str, indent: int) -> str:
    lines = doc.split("\n")
    for i, l in enumerate(lines):
        if len(l) == 0:
            continue
        lines[i] = (" " * indent) + l
    return "\n".join(lines)


def _members_to_doc(members: Sequence[MemberWithDoc]) -> str:
    content = ""
    for m in members:
        if m.name.startswith("_"):
            continue
        doc = m.docstring
        offset = 0
        sorted_directives = m.directives
        sorted_directives.sort(key=lambda k: k.start)
        for d in sorted_directives:
            doc_length = len(doc)
            if d.directive_type == LawDirectiveTypes.SYMBOL:
                directive_str = ":code:`" + code_str(m.value) + "`"
                doc = doc[0:d.start + offset] + directive_str + "\n" + doc[d.end + offset:]
                offset = offset + len(doc) - doc_length
                continue
            if d.directive_type == LawDirectiveTypes.LATEX:
                directive_str = "Latex:\n" + _indent_docstring(".. math::\n", INDENT_SPACES)
                directive_str = directive_str + _indent_docstring(latex_str(m.value), INDENT_SPACES * 2)
                doc = doc[0:d.start + offset] + directive_str + "\n" + doc[d.end + offset:]
                offset = offset + len(doc) - doc_length
                continue
        doc = _indent_docstring(doc, INDENT_SPACES)
        content = content + ".. py:data:: " + m.name + "\n\n"
        content = content + doc + "\n\n"
        if m.symbol is not None:
            symbol_name = "Symbol:\n" + _indent_docstring(":code:`" + m.symbol.symbol + "`", INDENT_SPACES)
            content = content + _indent_docstring(symbol_name, INDENT_SPACES) + "\n\n"
            if m.symbol.latex is not None:
                latex_string = m.symbol.latex
                symbol_latex = "Latex:\n" + _indent_docstring(":math:`" + latex_string + "`", INDENT_SPACES)
                content = content + _indent_docstring(symbol_latex, INDENT_SPACES) + "\n\n"
            symbol_dimension = "Dimension:\n" + _indent_docstring(":code:`" + m.symbol.dimension + "`", INDENT_SPACES)
            content = content + _indent_docstring(symbol_dimension, INDENT_SPACES) + "\n\n"
    return content


def _functions_to_doc(members: Sequence[FunctionWithDoc]) -> str:
    content = ""
    for m in members:
        if m.name.startswith("_"):
            continue
        doc = _indent_docstring(m.docstring, INDENT_SPACES)
        args = ", ".join(m.parameters)
        content = content + ".. py:function:: " + m.name + "(" + args + ")" + "\n\n"
        content = content + doc + "\n\n"
    return content


def print_law(title: str, description: str, items: Sequence[MemberWithDoc | FunctionWithDoc], doc_name: str) -> str:
    law_content = "" + title + "\n" + ("-" * len(title)) + "\n\n" + description + "\n\n"
    law_content = law_content + ".. py:currentmodule:: " + doc_name + "\n\n"
    members = [m for m in items if isinstance(m, MemberWithDoc)]
    functions = [m for m in items if isinstance(m, FunctionWithDoc)]
    law_content = law_content + _members_to_doc(members)
    law_content = law_content + _functions_to_doc(functions)

    return law_content


def print_package(title: str, description: str, items: Sequence[MemberWithDoc | FunctionWithDoc], doc_name: str, laws: Sequence[str], packages: Sequence[str]) -> str:
    package_content = "" + title + "\n" + ("=" * len(title)) + "\n\n" + description + "\n\n"
    package_content = package_content + ".. py:currentmodule:: " + doc_name + "\n\n"
    if len(packages) > 0 or len(laws) > 0:
        package_content = package_content + "Contents:\n\n" + ".. toctree::\n" + "  :maxdepth: 4\n\n"
        for p in packages:
            package_content = package_content + "  " + p + "\n"
        for l in laws:
            package_content = package_content + "  " + l + "\n"

    members = [m for m in items if isinstance(m, MemberWithDoc)]
    functions = [m for m in items if isinstance(m, FunctionWithDoc)]
    package_content = package_content + _members_to_doc(members)
    package_content = package_content + _functions_to_doc(functions)

    return package_content
