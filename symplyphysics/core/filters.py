# mypy: disable-error-code="arg-type, comparison-overlap"

from typing import Any, Mapping, Sequence, TypeVar
from abc import abstractmethod

K = TypeVar("K")
T = TypeVar("T", bound="Comparable")


class Comparable:

    @abstractmethod
    def __eq__(self, other: Any) -> bool:
        pass

    @abstractmethod
    def __lt__(self: T, other: T | int) -> bool:
        pass

    def __gt__(self: T, other: T | int) -> bool:
        return (not self < other) and self != other

    def __le__(self: T, other: T | int) -> bool:
        return self < other or self == other

    def __ge__(self: T, other: T | int) -> bool:
        return not self < other


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
