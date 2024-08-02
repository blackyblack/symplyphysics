from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.definitions import number_density_is_number_of_objects_per_unit_volume as law

Args = namedtuple("Args", ["o", "V"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    objects = 100
    V = Quantity(1 * units.meter**3)
    return Args(o=objects, V=V)


def test_basic_density(test_args: Args) -> None:
    result = law.calculate_number_density(test_args.o, test_args.V)
    assert_equal(result, 100 / units.meter**3)


def test_bad_volume(test_args: Args) -> None:
    Vb = Quantity(1 * units.length)
    with raises(errors.UnitsError):
        law.calculate_number_density(test_args.o, Vb)
    with raises(TypeError):
        law.calculate_number_density(test_args.o, 100)
