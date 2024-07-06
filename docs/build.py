#!/usr/bin/env python3

import argparse
import ast
import importlib
import os
import sys
from typing import Optional, Sequence
from sphinx.pycode import ModuleAnalyzer
from sphinx.application import Sphinx


def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        usage='%(prog)s [OPTIONS] <SOURCE_FILE>...',
        description='Generate ReStructuredText for Symplyphysics laws.'
    )

    parser.add_argument(
        '-s',
        '--source-dir',
        action='store',
        dest='source_dir',
        help='source directory to generate rST files for'
    )

    parser.add_argument(
        '-o',
        '--output-dir',
        action='store',
        dest='output_dir',
        help='directory to place all output in',
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


def process_law_package(directory: str, laws: Sequence[str], output_dir: str) -> Optional[str]:
    _, directory_name = os.path.split(directory)
    if directory_name.startswith("__"):
        return None
    filename = os.path.normpath(directory)
    filename_init = os.path.join(filename, "__init__.py")
    init_content: Optional[str] = None
    with open(filename_init, "r", encoding="utf-8") as f:
        init_content = f.read()
    if init_content is None:
        return None
    docstring = ast.get_docstring(ast.parse(init_content))
    if docstring is None:
        return None
    law_title = find_title(docstring)
    if law_title is None:
        return None
    law_description = find_description(docstring)
    law_description = "" if law_description is None else law_description

    package_doc_name = filename.replace(os.sep, ".") + ".rst"
    package_doc_file = os.path.normpath(os.path.join(output_dir, package_doc_name))

    package_content = "" + law_title + "\n" + ("=" * len(law_title)) + "\n\n" + law_description + "\n\n"
    package_content = package_content + "Contents:\n\n" + ".. toctree::\n" + "  :maxdepth: 4\n\n"
    for l in laws:
        package_content = package_content + "  " + l + "\n"

    # create directory if necessary
    os.makedirs(os.path.dirname(package_doc_file), exist_ok=True)
    with open(package_doc_file, "w+", encoding="utf-8") as f:
        f.write(package_content)
    return package_doc_name


def process_law(directory: str, law_filename: str, output_dir: str) -> Optional[str]:
    if not law_filename.endswith(".py"):
        return None
    if law_filename.startswith("__"):
        return None
    filename = os.path.normpath(os.path.join(directory, law_filename))
    no_extension_filename = os.path.splitext(filename)[0]
    law_module_name = no_extension_filename.replace(os.sep, ".")
    # TODO: this is very slow
    #module = importlib.import_module(law_module_name)
    #docstring = module.__doc__
    #if docstring is None:
    #    return None
    #law_title = find_title(docstring)
    #if law_title is None:
    #    return None
    #law_description = find_description(docstring)
    #law_description = "" if law_description is None else law_description
    law_title = "x"
    law_description = "y"

    # TODO: this is not very fast as well
    analyzer = ModuleAnalyzer.for_module(law_module_name)
    law_members: list[str] = []
    for _ns, name in analyzer.find_attr_docs():
        law_members.append(name)

    law_doc_name = law_module_name + ".rst"
    law_doc_file = os.path.normpath(os.path.join(output_dir, law_doc_name))

    law_content = "" + law_title + "\n" + ("-" * len(law_title)) + "\n\n" + law_description + "\n\n"
    for m in law_members:
        law_content = law_content + ".. autodata:: " + law_module_name + "." + m + "\n"
        law_content = law_content + "  :no-value:\n\n"

    # create directory if necessary
    os.makedirs(os.path.dirname(law_doc_file), exist_ok=True)
    with open(law_doc_file, "w+", encoding="utf-8") as f:
        f.write(law_content)
    return law_module_name


def generate_laws_docs(source_dir: str, output_dir: str) -> None:
    for path, _dirs, files in os.walk(source_dir):
        # TODO: underlying packages should also be included in packages

        laws: list[str] = []
        for filename in files:
            law_name = process_law(path, filename, output_dir)
            if law_name is None:
                continue
            laws.append(law_name)
        package_name = process_law_package(path, laws, output_dir)
        if package_name is None:
            continue


def main(argv: Sequence[str] = (), /) -> None:
    args = get_parser().parse_args(argv or sys.argv[1:])

    generate_laws_docs(args.source_dir, args.output_dir)

    # Build HTML docs
    app = Sphinx("docs", "docs", "html", "html/.doctrees", "html")
    app.build()


if __name__ == '__main__':
    main(sys.argv[1:])
