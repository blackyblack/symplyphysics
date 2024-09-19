from pathlib import Path
import re

from symplyphysics import quantities, Quantity

_role_pattern = re.compile(r":quantity_notation:`(\w*)`")

_collection = set()
for _attr in dir(quantities):
    _obj = getattr(quantities, _attr)
    if isinstance(_obj, Quantity):
        _collection.add(_attr)

def process_string(doc: str, path: Path) -> str:
    parts = []

    last_index_to = 0
    for match in _role_pattern.finditer(doc):
        index_from, index_to = match.span()
        name = match.group(1)

        if name not in _collection:
            raise ValueError(f"Unknown quantity {name} in '{path}'.")

        part_before = doc[last_index_to:index_from]

        qty: Quantity = getattr(quantities, name)
        text = f":math:`{qty.display_latex}` (:code:`{qty.display_name}`) is :attr:`~symplyphysics.quantities.{name}`"

        parts.append(part_before)
        parts.append(text)

        last_index_to = index_to

    if not parts:
        return doc

    part_after = doc[last_index_to:]
    parts.append(part_after)

    return "".join(parts)
