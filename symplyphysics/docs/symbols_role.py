from pathlib import Path
from types import ModuleType

import re

from symplyphysics import symbols, SymbolNew

_symbols_pattern = re.compile(r":symbols:`(\w*)`")

# Automatically find the necessaty modules so that we don't need to write the logic by hand
_symbols_by_module: dict[str, set[str]] = {}
for _attr in set(dir(symbols)) - set(symbols.__all__):
    _obj = getattr(symbols, _attr)
    _file = getattr(_obj, "__file__", "")
    if (
        isinstance(_obj, ModuleType)
        and Path(symbols.__file__).parent == Path(_file).parent
    ):
        _symbols_by_module[_attr] = set()
        for _subattr in dir(_obj):
            _subobj = getattr(_obj, _subattr)
            if isinstance(_subobj, SymbolNew):
                _symbols_by_module[_attr].add(_subattr)


def process_string(doc: str, path: str) -> str:
    parts = []

    last_index_to = 0
    for match in _symbols_pattern.finditer(doc):
        index_from, index_to = match.span()
        name = match.group(1)

        for directory, collection in _symbols_by_module.items():
            if name in collection:
                part_before = doc[last_index_to:index_from]
                text = f":attr:`~symplyphysics.symbols.{directory}.{name}`"

                parts.append(part_before)
                parts.append(text)

                last_index_to = index_to
                break
        else:
            raise ValueError(f"Unknown symbol {name} in '{path}'.")

    # no substitution happened
    if not parts:
        return doc

    part_after = doc[last_index_to:]
    parts.append(part_after)

    return "".join(parts)
