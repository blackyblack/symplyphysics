from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, units, Quantity, errors, prefixes, quantities)
from symplyphysics.laws.electricity import magnetic_field_due_to_infinite_wire as law

Args = namedtuple("Args", ["permeability", "current", "distance"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    permeability = Quantity(quantities.vacuum_permeability * 1.1)
    current = Quantity(0.5 * units.ampere)
    distance = Quantity(2 * units.meter)
    return Args(permeability=permeability, current=current, distance=distance)


def test_law(test_args: Args) -> None:
    result = law.calculate_magnetic_field(test_args.permeability, test_args.current,
        test_args.distance)
    assert_equal(result, 55 * prefixes.nano * units.tesla)


def test_bad_permeability(test_args: Args) -> None:
    mub = Quantity(units.candela)
    with raises(errors.UnitsError):
        law.calculate_magnetic_field(mub, test_args.current, test_args.distance)
    with raises(TypeError):
        law.calculate_magnetic_field(100, test_args.current, test_args.distance)


def test_bad_current(test_args: Args) -> None:
    current = Quantity(1 * units.candela)
    with raises(errors.UnitsError):
        law.calculate_magnetic_field(test_args.permeability, current, test_args.distance)
    with raises(TypeError):
        law.calculate_magnetic_field(test_args.permeability, 100, test_args.distance)


def test_bad_distance(test_args: Args) -> None:
    distance = Quantity(1 * units.candela)
    with raises(errors.UnitsError):
        law.calculate_magnetic_field(test_args.permeability, test_args.current, distance)
    with raises(TypeError):
        law.calculate_magnetic_field(test_args.permeability, test_args.current, 100)
