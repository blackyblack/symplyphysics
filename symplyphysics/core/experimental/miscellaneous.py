from typing import TypeVar, Iterable, Callable, Any, Optional
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
