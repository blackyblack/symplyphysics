from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, units, Quantity, errors, prefixes)
from symplyphysics.laws.electricity import magnetic_field_due_to_infinite_wire as law

Args = namedtuple("Args", ["current", "distance"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    current = Quantity(0.5 * units.ampere)
    distance = Quantity(2 * units.meter)
    return Args(current=current, distance=distance)


def test_law(test_args: Args) -> None:
    result = law.calculate_magnetic_field(test_args.current, test_args.distance)
    assert_equal(result, 50 * prefixes.nano * units.tesla)


def test_bad_current(test_args: Args) -> None:
    current = Quantity(1 * units.candela)
    with raises(errors.UnitsError):
        law.calculate_magnetic_field(current, test_args.distance)
    with raises(TypeError):
        law.calculate_magnetic_field(100, test_args.distance)


def test_bad_distance(test_args: Args) -> None:
    distance = Quantity(1 * units.candela)
    with raises(errors.UnitsError):
        law.calculate_magnetic_field(test_args.current, distance)
    with raises(TypeError):
        law.calculate_magnetic_field(test_args.current, 100)
