from __future__ import annotations

from typing import Optional
from pathlib import Path
from dataclasses import dataclass
from collections import defaultdict

import re

_HEADING_PAT = re.compile(r"^([^\n]+)")

_LATEX_PAT = re.compile(
    # pylint: disable-next=line-too-long
    r"\.\. py:data:: (?:law|definition|condition)\s+:code:[^\n]+\s+Latex:\s+\.\. math::\n((?: {12}(?:[^\n]+)\n)+)"
)


@dataclass(slots=True, frozen=True)
class _Info:
    parents: str
    name: str
    heading: str
    latex: str

    @staticmethod
    def try_parse(path: Path) -> Optional[_Info]:
        with open(path, "r", encoding="utf-8") as file:
            source = file.read()

        m = _HEADING_PAT.search(source)
        if not m:
            return None
        heading = m.group(1)

        m = _LATEX_PAT.search(source)
        if not m:
            return None
        latex = m.group(1).strip()

        *parents, name = path.stem.split(".")

        return _Info(parents=".".join(parents), name=name, heading=heading, latex=latex)


_HEADING_TEMPLATE = """\
Formulas
========

{body}
"""

_DIRECTORY_TEMPLATE = """\
- :doc:`(doc) <{parents}>` :code:`{parents}`

    .. toggle::

        {body}"""

_MODULE_TEMPLATE = """\
- :doc:`(doc) <{parents}.{name}>` {heading}

    .. math::

        {latex}"""


def _print_tree(tree: dict[str, list[_Info]]) -> str:
    directories = []

    for parents, infos in tree.items():
        modules = []

        for info in infos:
            module = _MODULE_TEMPLATE.format(
                parents=parents,
                name=info.name,
                heading=info.heading,
                latex=info.latex,
            )

            modules.append(module)

        body = "\n".join(modules).replace("\n", "\n        ")

        directory = _DIRECTORY_TEMPLATE.format(
            parents=parents,
            body=body,
        )

        directories.append(directory)

    return "\n".join(directories)


def generate_formulas(generated_dir: str, output_name: str) -> None:
    generated_path = Path(generated_dir)

    tree: dict[str, list[_Info]] = defaultdict(list)

    for path in generated_path.glob("*.rst"):
        info = _Info.try_parse(path)
        if not info:
            continue

        tree[info.parents].append(info)

    for parents in tree:
        tree[parents].sort(key=lambda info: info.name)

    tree = dict(sorted(tree.items(), key=lambda p: p[0]))

    body = _print_tree(tree)
    source = _HEADING_TEMPLATE.format(body=body)

    with open(generated_path / output_name, "w", encoding="utf-8") as file:
        file.write(source)
