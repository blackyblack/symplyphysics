from typing import TypeVar, Iterable, Callable, Any, Optional, Generic
from sympy.combinatorics.permutations import Permutation

_T = TypeVar("_T")


def sort_with_sign(
    it: Iterable[_T],
    key: Optional[Callable[[_T], Any]] = None,
) -> tuple[int, list[_T]]:
    """
    Sorts `it` with a optional `key`. Also returns the sign of the permutation of the elements of
    `it` as a result of sorting.
    """

    old_it = it if isinstance(it, list) else list(it)
    new_it = sorted(old_it, key=key)

    indeces = [old_it.index(v) for v in new_it]

    if len(set(indeces)) != len(indeces):
        sign = 0
    else:
        sign = Permutation(indeces).signature()

    return sign, new_it


class Registry(Generic[_T]):
    """
    A registry contains a mapping of the given elements to integers in the order of their addition
    into the registry.
    """

    _index: int
    _mapping: dict[_T, int]

    def __init__(self) -> None:
        self._index = 0
        self._mapping = {}

    def add(self, value: _T) -> None:
        self._index += 1
        self._mapping[value] = self._index

    def get(self, value: _T) -> int:
        return self._mapping[value]

    def __len__(self) -> int:
        return self._index
