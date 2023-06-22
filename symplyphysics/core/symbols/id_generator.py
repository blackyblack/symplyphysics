# Mapping from base prefix to the last assigned id
_ids: dict[str, int] = {}

# Assign and get next id
def next_id(base: str = "") -> int:
    id = _ids.get(base)
    id = 1 if id is None else id + 1
    _ids[base] = id
    return id
