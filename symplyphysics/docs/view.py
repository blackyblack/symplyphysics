"""
This module provides a functionality for printing laws and packages.

* `print_law` constructs the documentation for a law.

* `print_package` constructs the documentation for a package.
"""

from typing import Sequence, Optional
from .printer_latex import latex_str
from .printer_code import code_str
from .parse import FunctionWithDoc, LawDirectiveType, MemberWithDoc

_INDENT = "    "


def _indent_docstring(doc: str, count: int = 1) -> str:
    lines = doc.splitlines()
    for i, line in enumerate(lines):
        if not line:
            continue
        lines[i] = (_INDENT * count) + line
    return "\n".join(lines)


_SYMBOL_DIRECTIVE_TEMPLATE = """\
{before}:code:`{code}`
{after}\
"""

_LATEX_DIRECTIVE_TEMPLATE = """\
{before}Latex:
    .. math::
{ii_latex}
{after}\
"""

_MEMBER_TEMPLATE = """\
.. py:data:: {name}

{i_doc}
"""

_LAW_SYMBOL_TEMPLATE = """\
Symbol:
    :code:`{code}`

Latex:
    :math:`{latex}`

Dimension:
    :code:`{dimension}`
"""


def _members_to_doc(members: Sequence[MemberWithDoc], doc_name: str) -> str:

    def process_member_docstring(member: MemberWithDoc) -> str:
        doc = member.docstring
        offset = 0

        for directive in sorted(member.directives, key=lambda k: k.start):
            doc_length = len(doc)

            before = doc[:directive.start + offset]
            after = doc[directive.end + offset:]

            match directive.directive_type:
                case LawDirectiveType.SYMBOL:
                    try:
                        code = code_str(member.value)
                    except ValueError as e:
                        raise ValueError(f"Error during code printing in '{doc_name}'.") from e

                    doc = _SYMBOL_DIRECTIVE_TEMPLATE.format(
                        before=before,
                        code=code,
                        after=after,
                    )

                case LawDirectiveType.LATEX:
                    try:
                        latex = latex_str(member.value)
                    except ValueError as e:
                        raise ValueError(f"Error during latex printing in '{doc_name}'.") from e

                    doc = _LATEX_DIRECTIVE_TEMPLATE.format(
                        before=before,
                        ii_latex=_indent_docstring(latex, count=2),
                        after=after,
                    )

            offset += len(doc) - doc_length

        return doc

    def member_to_doc(member: MemberWithDoc) -> Optional[str]:
        if member.name.startswith("_"):
            return None

        docstring = process_member_docstring(member)
        doc = _MEMBER_TEMPLATE.format(name=member.name, i_doc=_indent_docstring(docstring))

        if member.symbol:
            symbol = _LAW_SYMBOL_TEMPLATE.format(
                code=member.symbol.symbol,
                latex=member.symbol.latex,
                dimension=member.symbol.dimension,
            )
            doc += f"\n{symbol}"

        return doc

    docs = []
    for member in members:
        doc = member_to_doc(member)
        if doc is not None:
            docs.append(doc)

    return "\n\n".join(docs)


_FUNCTION_TEMPLATE = """\
.. py:function:: {name}({args})

{i_doc}
"""


def _functions_to_doc(members: Sequence[FunctionWithDoc]) -> str:

    def function_to_doc(function_: FunctionWithDoc) -> Optional[str]:
        if function_.name.startswith("_"):
            return None

        return _FUNCTION_TEMPLATE.format(
            name=function_.name,
            args=", ".join(function_.parameters),
            i_doc=_indent_docstring(function_.docstring),
        )

    docs = []
    for function_ in members:
        doc = function_to_doc(function_)
        if doc is not None:
            docs.append(doc)

    return "\n\n".join(docs)


_HEADER_TEMPLATE = """\
{title}
{section_break}

{description}

.. py:currentmodule:: {doc_name}\
"""

_CONTENTS_TEMPLATE = """\
**Contents:**

.. toctree::
    :maxdepth: 4

{i_packages}
{i_laws}\
"""

_FOOTER_TEMPLATE = """\
{members}

{functions}\
"""


def print_law(title: str, description: str, members: Sequence[MemberWithDoc],
    functions: Sequence[FunctionWithDoc], doc_name: str) -> str:
    """Constructs the documentation string of a law."""

    header = _HEADER_TEMPLATE.format(
        title=title,
        section_break="=" * len(title),
        description=description,
        doc_name=doc_name,
    )

    footer = _FOOTER_TEMPLATE.format(
        members=_members_to_doc(members, doc_name),
        functions=_functions_to_doc(functions),
    )

    return f"{header}\n\n{footer}"


# TODO: split to smaller functions


# pylint: disable-next=too-many-arguments, too-many-positional-arguments
def print_package(title: str, description: str, members: Sequence[MemberWithDoc],
    functions: Sequence[FunctionWithDoc], doc_name: str, laws: Sequence[str],
    packages: Sequence[str]) -> str:
    """Construct the documentation string of a package."""

    results = []

    header = _HEADER_TEMPLATE.format(
        title=title,
        section_break="=" * len(title),
        description=description,
        doc_name=doc_name,
    )
    results.append(header)

    if packages or laws:
        joined_packages = "\n".join(map(_indent_docstring, packages))
        joined_laws = "\n".join(map(_indent_docstring, laws))
        contents = _CONTENTS_TEMPLATE.format(
            i_packages=joined_packages,
            i_laws=joined_laws,
        )
        results.append(contents)

    footer = _FOOTER_TEMPLATE.format(
        members=_members_to_doc(members, doc_name),
        functions=_functions_to_doc(functions),
    )
    results.append(footer)

    return "\n\n".join(results)


__all__ = ["print_law", "print_package"]
