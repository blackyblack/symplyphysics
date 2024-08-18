import ast
import os
from typing import Optional, Sequence
from .parse import find_description, find_members_and_functions, find_title
from .view import print_law, print_package


def process_law_package(directory: str, laws: Sequence[str], packages: Sequence[str],
    output_dir: str) -> Optional[str]:
    filename = os.path.normpath(directory)
    filename_init = os.path.join(filename, "__init__.py")
    package_content: Optional[str] = None
    with open(filename_init, "r", encoding="utf-8") as f:
        package_content = f.read()
    package_parsed = ast.parse(package_content)
    docstring = ast.get_docstring(package_parsed)
    if docstring is None:
        return None
    package_title = find_title(docstring)
    if package_title is None:
        return None
    package_description = find_description(docstring)
    if package_description is None:
        package_description = ""

    package_name = filename.replace(os.sep, ".")
    package_name_strip_root = ".".join(package_name.split(".")[1:])
    packages = [p for p in packages if not (p.startswith(".") or p.startswith("_"))]
    packages = [package_name_strip_root + "." + p for p in packages]
    package_items = find_members_and_functions(package_parsed)
    package_content = print_package(package_title, package_description, package_items, package_name, laws, packages)

    package_doc_file = os.path.normpath(os.path.join(output_dir, package_name_strip_root + ".rst"))

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
    if law_description is None:
        law_description = ""

    no_extension_filename = os.path.splitext(filename)[0]
    law_module_name = no_extension_filename.replace(os.sep, ".")
    law_items = find_members_and_functions(law_parsed)
    doc_content = print_law(law_title, law_description, law_items, law_module_name)

    law_module_name_strip_root = ".".join(law_module_name.split(".")[1:])
    doc_name = law_module_name_strip_root + ".rst"
    doc_file = os.path.normpath(os.path.join(output_dir, doc_name))

    # create directory if necessary
    os.makedirs(os.path.dirname(doc_file), exist_ok=True)
    with open(doc_file, "w+", encoding="utf-8") as f:
        f.write(doc_content)
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
        dirs.sort()
        files.sort()
        laws: list[str] = []
        for filename in files:
            law_name = process_law(path, filename, output_dir)
            if law_name is None:
                continue
            laws.append(law_name)
        package_name = process_law_package(path, laws, dirs, output_dir)
        if package_name is None:
            continue
