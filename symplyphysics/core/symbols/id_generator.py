"""
This module maintains a global list of generated IDs.
"""

# Mapping from base prefix to the last assigned id
_ids: dict[str, int] = {}


# Assign and get next id
def next_id(base: str = "") -> int:
    id_val = _ids.get(base)
    id_val = 1 if id_val is None else id_val + 1
    _ids[base] = id_val
    return id_val
