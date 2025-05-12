from typing import TypeVar, Iterable, Callable, Any, Optional, Sequence, Generator
from sympy import cacheit as sym_cacheit, Expr, sympify
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


def cacheit(f: _T) -> _T:
    """A typed version of `sympy.cacheit`."""

    return sym_cacheit(f)


def select_by_indices(items: Sequence[_T], indices: Iterable[int]) -> Generator[_T, None, None]:
    for index in indices:
        yield items[index]


def sympify_expr(value: Any) -> Expr:
    result = sympify(value, strict=True)

    if not isinstance(result, Expr):
        raise TypeError(f"Expected '{result}' to be Expr.")

    return result
