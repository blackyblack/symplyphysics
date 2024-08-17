import ast
from dataclasses import dataclass
from typing import Optional


@dataclass
class MemberWithDoc:
    name: str
    docstring: str

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


def find_members_and_functions(content: ast.Module) -> list[MemberWithDoc | FunctionWithDoc]:
    law_functions: list[FunctionWithDoc] = []
    law_members: list[str] = []
    current_member: Optional[str] = None
    docstrings: dict[str, str] = {}
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
        result.append(MemberWithDoc(v, doc))
    for f in law_functions:
        result.append(f)
    return result
