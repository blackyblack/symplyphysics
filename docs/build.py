#!/usr/bin/env python3

import argparse
import shutil
import sys
from typing import Sequence
from pathlib import Path
from sphinx.application import Sphinx

from symplyphysics.docs.build import generate_laws_docs
from symplyphysics.docs import symbols_role


def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(usage="%(prog)s [OPTIONS] <SOURCE_FILE>...",
        description="Generate ReStructuredText for Symplyphysics laws.")

    parser.add_argument("-l",
        "--laws-source-dir",
        action="store",
        dest="laws_source_dir",
        default="symplyphysics",
        help="laws source directory to generate rST files from")

    parser.add_argument(
        "-g",
        "--generated-dir",
        action="store",
        dest="generated_dir",
        default="docs/generated",
        help="directory to place all generated files",
    )

    parser.add_argument(
        "-e",
        "--exclude-dirs",
        nargs="*",
        default=["core"],
        help="directories to exclude from parsing",
    )

    parser.add_argument("-o",
        "--output-dir",
        action="store",
        dest="output_dir",
        default="html",
        help="directory to generate HTML output into")

    parser.add_argument("-c",
        "--conf-dir",
        action="store",
        dest="conf_dir",
        default="docs",
        help="directory where conf.py file is stored")
    
    parser.add_argument("-q",
        "--quiet",
        action="store_true",
        dest="quiet",
        help="suppress rST files generation text",
        )

    return parser


def process_generated_files(generated_dir: str) -> None:
    for file_path in Path(generated_dir).iterdir():
        with open(file_path, "r+", encoding="utf-8") as file:
            doc = file.read()
            
            # processing logic goes here
            doc = symbols_role.process_string(doc)

            file.write(doc)


def main(argv: Sequence[str] = ()) -> None:
    args = get_parser().parse_args(argv or sys.argv[1:])

    # Generate Symplyphysics rst files
    generate_laws_docs(args.laws_source_dir, args.generated_dir, args.exclude_dirs, args.quiet)

    # Copy index.rst to 'generated' folder
    index_file = Path(args.conf_dir) / "index.rst"
    out_index_file = Path(args.generated_dir) / "index.rst"
    shutil.copyfile(index_file, out_index_file, follow_symlinks=True)

    process_generated_files(args.generated_dir)

    # Build HTML docs
    doctrees = Path(args.output_dir) / ".doctrees/"
    app = Sphinx(args.generated_dir, args.conf_dir, "html", doctrees, args.output_dir)
    app.build()


if __name__ == "__main__":
    main(sys.argv[1:])
