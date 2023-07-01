from typing import Any, Mapping, Sequence, TypeVar


class Comparable:

    def __gt__(self, _: Any) -> bool:
        ...

    def __lt__(self, _: Any) -> bool:
        ...

    def __ge__(self, _: Any) -> bool:
        ...

    def __le__(self, _: Any) -> bool:
        ...

    def __eq__(self, _: Any) -> bool:
        ...


T = TypeVar("T", bound=Comparable)
K = TypeVar("K")


def filter_zeroes(list_: Sequence[T]) -> list[T]:
    return list(filter(lambda e: e != 0, list_))


# list_ contains map as an element
def filter_map_zeroes(key_: K, list_: Sequence[Mapping[K, T]]) -> list[Mapping[K, T]]:
    return list(filter(lambda e: e[key_] != 0, list_))


def filter_negative(list_: Sequence[T]) -> list[T]:
    return list(filter(lambda e: e >= 0, list_))


# list_ contains map as an element
def filter_map_negative(key_: K, list_: Sequence[Mapping[K, T]]) -> list[Mapping[K, T]]:
    return list(filter(lambda e: e[key_] >= 0, list_))
