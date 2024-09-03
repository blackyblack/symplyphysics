#!/usr/bin/env python3

import argparse
import pathlib
import shutil
import sys
from typing import Sequence
from sphinx.application import Sphinx

from symplyphysics.docs.build import generate_laws_docs


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

    return parser


def main(argv: Sequence[str] = ()) -> None:
    args = get_parser().parse_args(argv or sys.argv[1:])

    # Generate Symplyphysics rst files
    generate_laws_docs(args.laws_source_dir, args.generated_dir, args.exclude_dirs)

    # Copy index.rst to 'generated' folder
    index_file = pathlib.Path(args.conf_dir).joinpath("index.rst")
    out_index_file = pathlib.Path(args.generated_dir).joinpath("index.rst")
    shutil.copyfile(index_file, out_index_file, follow_symlinks=True)

    # Build HTML docs
    doctrees = pathlib.Path(args.output_dir).joinpath(".doctrees/")
    app = Sphinx(args.generated_dir, args.conf_dir, "html", doctrees, args.output_dir)
    app.build()


if __name__ == "__main__":
    main(sys.argv[1:])
