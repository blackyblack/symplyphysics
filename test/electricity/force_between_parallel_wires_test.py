from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, units, Quantity, errors, quantities)
from symplyphysics.laws.electricity import force_between_parallel_wires as force_law

# Description
## Let the current in the first wire be 0.5 amperes, and in the second 1.5 amperes.
## Then, with a length of wires equal to 1 meter and a relative permeability equal to 1,
## the force will be equal to 7.5e-8 newtons at a distance of 2 meters.
## https://www.calculatoratoz.com/ru/force-between-parallel-wires-calculator/Calc-2131

Args = namedtuple("Args",
    ["absolute_permeability", "first_wire_current", "second_wire_current", "length", "distance"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    relative_permeability_ = 1
    absolute_permeability = Quantity(relative_permeability_ * quantities.vacuum_permeability)
    first_wire_current = Quantity(0.5 * units.ampere)
    second_wire_current = Quantity(1.5 * units.ampere)
    length = Quantity(1 * units.meter)
    distance = Quantity(2 * units.meter)
    return Args(absolute_permeability=absolute_permeability,
        first_wire_current=first_wire_current,
        second_wire_current=second_wire_current,
        length=length,
        distance=distance)


def test_basic_force(test_args: Args) -> None:
    result = force_law.calculate_force(test_args.absolute_permeability,
        test_args.first_wire_current, test_args.second_wire_current, test_args.length,
        test_args.distance)
    assert_equal(result, 7.5e-8 * units.newton)


def test_bad_absolute_permeability(test_args: Args) -> None:
    bad_absolute_permeability = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        force_law.calculate_force(bad_absolute_permeability, test_args.first_wire_current,
            test_args.second_wire_current, test_args.length, test_args.distance)
    with raises(TypeError):
        force_law.calculate_force(100, test_args.first_wire_current,
            test_args.second_wire_current, test_args.length, test_args.distance)


def test_bad_wire_current(test_args: Args) -> None:
    first_wire_current = Quantity(1 * units.coulomb)
    second_wire_current = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        force_law.calculate_force(test_args.absolute_permeability, first_wire_current,
            test_args.second_wire_current, test_args.length, test_args.distance)
    with raises(TypeError):
        force_law.calculate_force(test_args.absolute_permeability, 100,
            test_args.second_wire_current, test_args.length, test_args.distance)
    with raises(errors.UnitsError):
        force_law.calculate_force(test_args.absolute_permeability, test_args.first_wire_current,
            second_wire_current, test_args.length, test_args.distance)
    with raises(TypeError):
        force_law.calculate_force(test_args.absolute_permeability, test_args.first_wire_current,
            100, test_args.length, test_args.distance)


def test_bad_length(test_args: Args) -> None:
    length = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        force_law.calculate_force(test_args.absolute_permeability, test_args.first_wire_current,
            test_args.second_wire_current, length, test_args.distance)
    with raises(TypeError):
        force_law.calculate_force(test_args.absolute_permeability, test_args.first_wire_current,
            test_args.second_wire_current, 100, test_args.distance)


def test_bad_distance(test_args: Args) -> None:
    distance = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        force_law.calculate_force(test_args.absolute_permeability, test_args.first_wire_current,
            test_args.second_wire_current, test_args.length, distance)
    with raises(TypeError):
        force_law.calculate_force(test_args.absolute_permeability, test_args.first_wire_current,
            test_args.second_wire_current, test_args.length, 100)
