from ..core.symbols.symbols import Symbol
from . import (
    basic,
    kinematic,
    dynamics,
    thermodynamics,
)

__all__ = ["basic", "kinematic", "dynamics", "thermodynamics"]

_all_symbols = set()
for k, v in list(globals().items()):
    if not k in __all__:
        continue
    if isinstance(v, str):
        continue
    for val_name, symbol in v.__dict__.items():
        if not isinstance(symbol, Symbol):
            continue
        assert val_name not in _all_symbols, f"Duplicate symbol '{val_name}' definition."
        _all_symbols.add(val_name)
