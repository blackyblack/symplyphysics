#!/usr/bin/env python3

import argparse
import ast
import os
import sys
from typing import Optional, Sequence
from sphinx.application import Sphinx


def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(usage='%(prog)s [OPTIONS] <SOURCE_FILE>...',
        description='Generate ReStructuredText for Symplyphysics laws.')

    parser.add_argument('-s',
        '--source-dir',
        action='store',
        dest='source_dir',
        help='source directory to generate rST files for')

    parser.add_argument(
        '-o',
        '--output-dir',
        action='store',
        dest='output_dir',
        help='directory to place all output in',
    )

    parser.add_argument(
        '-e',
        '--exclude-dirs',
        nargs="*",
        help='directory to exclude from parsing',
    )

    return parser


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


def find_members_and_functions(content: ast.Module) -> tuple[list[str], list[str]]:
    law_functions: list[str] = []
    law_members: list[str] = []
    current_member: Optional[str] = None
    docstrings: dict[str, str] = {}
    for e in content.body:
        if isinstance(e, ast.FunctionDef):
            law_functions.append(e.name)
            current_member = e.name
            doc = ast.get_docstring(e)
            if doc is not None:
                docstrings[current_member] = doc
            continue
        if isinstance(e, ast.Assign):
            for t in e.targets:
                name = getattr(t, "id")
                law_members.append(name)
                current_member = name
            continue
        if isinstance(e, ast.Expr) and current_member is not None and isinstance(
                e.value, ast.Constant):
            docstrings[current_member] = e.value.value
            continue
    return ([m for m in law_members if m in docstrings],
        [m for m in law_functions if m in docstrings])


def members_to_doc(members: Sequence[str], doc_name: str) -> str:
    content = ""
    for m in members:
        if m.startswith("_"):
            continue
        content = content + ".. autodata:: " + doc_name + "." + m + "\n"
        content = content + "  :no-value:\n\n"
    return content


def functions_to_doc(members: Sequence[str], doc_name: str) -> str:
    content = ""
    for m in members:
        if m.startswith("_"):
            continue
        content = content + ".. autofunction:: " + doc_name + "." + m + "\n\n"
    return content


def process_law_package(directory: str, laws: Sequence[str], packages: Sequence[str],
    output_dir: str) -> Optional[str]:
    filename = os.path.normpath(directory)
    filename_init = os.path.join(filename, "__init__.py")
    init_content: Optional[str] = None
    with open(filename_init, "r", encoding="utf-8") as f:
        init_content = f.read()
    law_parsed = ast.parse(init_content)
    docstring = ast.get_docstring(law_parsed)
    if docstring is None:
        return None
    law_title = find_title(docstring)
    if law_title is None:
        return None
    law_description = find_description(docstring)
    law_description = "" if law_description is None else law_description

    law_members, law_functions = find_members_and_functions(law_parsed)

    package_name = filename.replace(os.sep, ".")
    package_name_strip_root = ".".join(package_name.split(".")[1:])
    package_doc_file = os.path.normpath(os.path.join(output_dir, package_name_strip_root + ".rst"))

    packages = [p for p in packages if not (p.startswith(".") or p.startswith("_"))]

    package_content = "" + law_title + "\n" + ("=" *
        len(law_title)) + "\n\n" + law_description + "\n\n"
    if len(packages) > 0 or len(laws) > 0:
        package_content = package_content + "Contents:\n\n" + ".. toctree::\n" + "  :maxdepth: 4\n\n"
        for p in packages:
            package_content = package_content + "  " + package_name_strip_root + "." + p + "\n"
        for l in laws:
            package_content = package_content + "  " + l + "\n"

    package_content = package_content + members_to_doc(law_members, package_name)
    package_content = package_content + functions_to_doc(law_functions, package_name)

    # create directory if necessary
    os.makedirs(os.path.dirname(package_doc_file), exist_ok=True)
    with open(package_doc_file, "w+", encoding="utf-8") as f:
        f.write(package_content)
    return package_name_strip_root


def process_law(directory: str, law_filename: str, output_dir: str) -> Optional[str]:
    if not law_filename.endswith(".py"):
        return None
    if law_filename.startswith("__"):
        return None
    filename = os.path.normpath(os.path.join(directory, law_filename))
    no_extension_filename = os.path.splitext(filename)[0]
    law_module_name = no_extension_filename.replace(os.sep, ".")
    law_module_name_strip_root = ".".join(law_module_name.split(".")[1:])

    with open(filename, "r", encoding="utf-8") as f:
        law_content = f.read()
    law_parsed = ast.parse(law_content)
    docstring = ast.get_docstring(law_parsed)
    if docstring is None:
        return None
    law_title = find_title(docstring)
    if law_title is None:
        return None
    law_description = find_description(docstring)
    law_description = "" if law_description is None else law_description

    law_members, law_functions = find_members_and_functions(law_parsed)

    law_doc_name = law_module_name_strip_root + ".rst"
    law_doc_file = os.path.normpath(os.path.join(output_dir, law_doc_name))

    law_content = "" + law_title + "\n" + ("-" * len(law_title)) + "\n\n" + law_description + "\n\n"
    law_content = law_content + members_to_doc(law_members, law_module_name)
    law_content = law_content + functions_to_doc(law_functions, law_module_name)

    # create directory if necessary
    os.makedirs(os.path.dirname(law_doc_file), exist_ok=True)
    with open(law_doc_file, "w+", encoding="utf-8") as f:
        f.write(law_content)
    return law_module_name_strip_root


def generate_laws_docs(source_dir: str, output_dir: str, exclude_dirs: Sequence[str]) -> None:
    exclude_dirs_with_source: list[str] = []
    for e in exclude_dirs:
        to_exclude = os.path.join(source_dir, e)
        exclude_dirs_with_source.append(to_exclude)
    for path, dirs, files in os.walk(source_dir):
        _, target_dir = os.path.split(path)
        if target_dir.startswith(".") or target_dir.startswith("_"):
            dirs.clear()
            files.clear()
            continue
        should_exclude = False
        for e in exclude_dirs_with_source:
            if os.path.samefile(path, e):
                should_exclude = True
                break
        if should_exclude:
            dirs.clear()
            files.clear()
            continue
        laws: list[str] = []
        for filename in files:
            law_name = process_law(path, filename, output_dir)
            if law_name is None:
                continue
            laws.append(law_name)
        package_name = process_law_package(path, laws, dirs, output_dir)
        if package_name is None:
            continue


def main(argv: Sequence[str] = (), /) -> None:
    args = get_parser().parse_args(argv or sys.argv[1:])

    generate_laws_docs(args.source_dir, args.output_dir, args.exclude_dirs)

    # Build HTML docs
    app = Sphinx("docs", "docs", "html", "html/.doctrees", "html")
    app.build()


if __name__ == '__main__':
    main(sys.argv[1:])
