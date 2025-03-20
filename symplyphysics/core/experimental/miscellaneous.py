from typing import TypeVar, Iterable, Callable, Any, Optional, Generic
from sympy import cacheit as sym_cacheit
from sympy.combinatorics.permutations import Permutation

_T = TypeVar("_T")


def sort_with_sign(
    it: Iterable[_T],
    key: Optional[Callable[[_T], Any]] = None,
) -> tuple[int, list[_T]]:
    """
    Sorts `it` with an optional `key`. Also returns the sign of the permutation of the elements of
    `it` as the result of sorting, which is `1` if the number of swaps performed between any two
    elements is even, `-1` if the number is odd, and `0` if `it` contains any equal elements.
    """

    old_it = it if isinstance(it, list) else list(it)

    if key is None:
        new_it = sorted(old_it)
        indices = [old_it.index(v) for v in new_it]
    else:
        old_keys = [key(v) for v in old_it]
        new_keys = sorted(old_keys)
        indices = [old_keys.index(k) for k in new_keys]
        new_it = [old_it[i] for i in indices]

    if len(set(indices)) != len(indices):
        sign = 0
    else:
        sign = Permutation(indices).signature()

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
        """Registers `value` in `self` unless it is already included in the registry."""

        if value in self._mapping:
            return

        self._index += 1
        self._mapping[value] = self._index

    def get(self, value: _T) -> int:
        """
        Gets the index of `value`. Does not check if the value is absent in the registry, which
        would raise a `KeyError`.
        """

        return self._mapping[value]

    def __len__(self) -> int:
        """Returns the number of elements registered."""

        return self._index


_T_co = TypeVar("_T_co", covariant=True)
_F = TypeVar("_F", bound=Callable[..., _T_co])


def cacheit(f: _F) -> _F:
    """A typed version of `sympy.cacheit`."""

    return sym_cacheit(f)  # type: ignore[no-any-return]
