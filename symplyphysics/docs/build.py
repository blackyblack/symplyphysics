"""
This module provides the function `generate_laws_docs` that generates rST files for laws and
packages in a given directory.
"""

import ast
import os
from pathlib import Path
from typing import Optional, Sequence
from .parse import find_members_and_functions, find_title_and_description
from .view import print_law, print_package
from .patch import patch_sympy_evaluate


def _parse_documentation(source: ast.Module) -> Optional[tuple[str, str]]:
    """
    Returns\\:
     
    * the title and description parsed from the ``source`` docstring.

    * `None` if there is no docstring or title.
    """

    docstring = ast.get_docstring(source)
    if docstring is None:
        return None

    return find_title_and_description(docstring)


def _is_private(s: str) -> bool:
    """Checks if the given string begins with `.` or `_`."""

    return s.startswith(".") or s.startswith("_")


def _process_law_package(directory: str, laws: Sequence[str], packages: Sequence[str],
    output_dir: str, quiet: bool) -> Optional[str]:
    """
    Processes the `__init__.py` file of the law package ``directory`` containing ``laws`` and
    sub-``packages``.

    Writes the resulting rST file in ``output_dir``.

    Suppresses generating messages if ``quiet`` is `True`.

    Returns\\:

    * the stem of the path to the written package documentation file.

    * `None` if the source file doesn't have any documentation.
    """

    source_dirpath = Path(directory)

    if not quiet:
        print(f"Generating rST for package {source_dirpath}.")

    source_init_path = source_dirpath / "__init__.py"
    with open(source_init_path, "r", encoding="utf-8", newline=None) as file:
        source_content = file.read()

    parsed_source = ast.parse(source_content)
    parsed_doc = _parse_documentation(parsed_source)
    if parsed_doc is None:
        return None
    title, description = parsed_doc

    parsed_source = patch_sympy_evaluate(parsed_source)

    try:
        members, functions = find_members_and_functions(parsed_source)
    except Exception as e:
        raise ValueError(f"Exception has been raised in {source_dirpath}.") from e

    source_file_stem = ".".join(source_dirpath.parts)

    doc_file_stem = ".".join(source_dirpath.parts[1:])
    packages = [doc_file_stem + "." + p for p in packages if not _is_private(p)]
    doc_content = print_package(title, description, members, functions, source_file_stem, laws,
        packages)

    doc_filepath = Path(output_dir, doc_file_stem + ".rst")

    # create directory if necessary
    doc_filepath.parent.mkdir(exist_ok=True)
    with open(doc_filepath, "w+", encoding="utf-8") as file:
        file.write(doc_content)

    return doc_file_stem


def _process_law(directory: str, filename: str, output_dir: str, quiet: bool) -> Optional[str]:
    """
    Processes the law in the given ``directory`` under the ``filename``.

    Writes the resulting rST file in ``output_dir``.

    Suppresses generating messages if ``quiet`` is `True`.

    Returns\\:

    * the stem of the path to the law documentation file.

    * `None` if the source file doesn't have any documentation.
    """

    if filename.startswith("__") or not filename.endswith(".py"):
        return None

    source_filepath = Path(directory, filename)

    if not quiet:
        print(f"Generating rST for {source_filepath}.")

    with open(source_filepath, "r", encoding="utf-8", newline=None) as file:
        source_content = file.read()

    parsed_source = ast.parse(source_content)
    parsed_doc = _parse_documentation(parsed_source)
    if parsed_doc is None:
        return None
    title, description = parsed_doc

    parsed_source = patch_sympy_evaluate(parsed_source)

    try:
        members, functions = find_members_and_functions(parsed_source)
    except Exception as e:
        raise ValueError(f"Exception has been raised in {source_filepath}.") from e

    source_filepath_strip_ext = source_filepath.with_suffix("")
    doc_long_stem = ".".join(source_filepath_strip_ext.parts)
    doc_content = print_law(title, description, members, functions, doc_long_stem)

    doc_short_stem = ".".join(source_filepath_strip_ext.parts[1:])
    doc_filepath = Path(output_dir, doc_short_stem + ".rst")

    # create directory if necessary
    doc_filepath.parent.mkdir(exist_ok=True)
    with open(doc_filepath, "w+", encoding="utf-8") as file:
        file.write(doc_content)

    return doc_short_stem


def generate_laws_docs(source_dir: str, output_dir: str, exclude_dirs: Sequence[str],
    quiet: bool) -> None:
    """
    Generates rST files for laws and packages, reading recursively starting from ``source_dir``,
    but avoiding ``exclude_dirs``.

    Writes the resulting files in ``output_dir``.
    """

    exclude_dirs_paths = [Path(source_dir, e) for e in exclude_dirs]

    for path, dirs, files in os.walk(source_dir):
        path = Path(path)

        target_dir = path.name
        if (_is_private(target_dir) or any(path.samefile(e) for e in exclude_dirs_paths)):
            # prevent `os.walk` from iterating excluded sub-directories
            dirs.clear()
            files.clear()
            continue

        dirs.sort()
        files.sort()

        laws: list[str] = []
        for file in files:
            law_name = _process_law(path, file, output_dir, quiet)
            if law_name is None:
                continue
            laws.append(law_name)

        package_name = _process_law_package(path, laws, dirs, output_dir, quiet)
        if package_name is None:
            continue


__all__ = ["generate_laws_docs"]
