from pathlib import Path
from types import ModuleType

import re

from symplyphysics import symbols, SymbolNew

_symbols_pattern = re.compile(r":symbols:`([^`]*)`")

_symbols_by_module: dict[str, set[str]] = {}
for _attr in set(dir(symbols)) - set(symbols.__all__):
    _obj = getattr(symbols, _attr)
    if (
        isinstance(_obj, ModuleType)
        and Path(symbols.__file__).parent == Path(_obj.__file__).parent
    ):
        _symbols_by_module[_attr] = set()
        for _subattr in dir(_obj):
            _subobj = getattr(_obj, _subattr)
            if isinstance(_subobj, SymbolNew):
                _symbols_by_module[_attr].add(_subattr)


def process_string(doc: str) -> str:
    matches = tuple(_symbols_pattern.finditer(doc))
    parts = []

    for match in reversed(matches):
        index_from, index_to = match.span()
        name = match.group(1)

        directory = None
        for the_directory, the_symbols in _symbols_by_module.items():
            if name in the_symbols:
                directory = the_directory
        if directory is None:
            raise ValueError(f"Unknown symbol {name}.")

        text = f":attr:`~symplyphysics.symbols.{directory}.{name}`"

        parts.append(doc[index_to:])
        parts.append(text)
        doc = doc[:index_from]

    return "".join(reversed(parts))
