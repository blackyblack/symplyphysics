from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, units, Quantity, errors)
from symplyphysics.laws.electricity import force_between_parallel_wires as force_law

# Description
## Let the current in the first wire be 0.5 amperes, and in the second 1.5 amperes.
## Then, with a length of wires equal to 1 meter and a relative permeability equal to 1,
## the force will be equal to 7.5e-8 newtons at a distance of 2 meters.
## https://www.calculatoratoz.com/ru/force-between-parallel-wires-calculator/Calc-2131


@fixture(name="test_args")
def test_args_fixture():
    relative_permeability = 1
    first_wire_current = Quantity(0.5 * units.ampere)
    second_wire_current = Quantity(1.5 * units.ampere)
    length = Quantity(1 * units.meter)
    distance = Quantity(2 * units.meter)

    Args = namedtuple("Args", [
        "relative_permeability", "first_wire_current", "second_wire_current", "length", "distance"
    ])
    return Args(relative_permeability=relative_permeability,
        first_wire_current=first_wire_current,
        second_wire_current=second_wire_current,
        length=length,
        distance=distance)


def test_basic_force(test_args):
    result = force_law.calculate_force(test_args.relative_permeability,
        test_args.first_wire_current, test_args.second_wire_current, test_args.length,
        test_args.distance)
    assert_equal(result, 7.5e-8 * units.newton)


def test_bad_relative_permeability(test_args):
    relative_permeability = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        force_law.calculate_force(relative_permeability, test_args.first_wire_current,
            test_args.second_wire_current, test_args.length, test_args.distance)


def test_bad_wire_current(test_args):
    first_wire_current = Quantity(1 * units.coulomb)
    second_wire_current = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        force_law.calculate_force(test_args.relative_permeability, first_wire_current,
            test_args.second_wire_current, test_args.length, test_args.distance)
    with raises(TypeError):
        force_law.calculate_force(test_args.relative_permeability, 100,
            test_args.second_wire_current, test_args.length, test_args.distance)
    with raises(errors.UnitsError):
        force_law.calculate_force(test_args.relative_permeability, test_args.first_wire_current,
            second_wire_current, test_args.length, test_args.distance)
    with raises(TypeError):
        force_law.calculate_force(test_args.relative_permeability, test_args.first_wire_current,
            100, test_args.length, test_args.distance)


def test_bad_length(test_args):
    length = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        force_law.calculate_force(test_args.relative_permeability, test_args.first_wire_current,
            test_args.second_wire_current, length, test_args.distance)
    with raises(TypeError):
        force_law.calculate_force(test_args.relative_permeability, test_args.first_wire_current,
            test_args.second_wire_current, 100, test_args.distance)


def test_bad_distance(test_args):
    distance = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        force_law.calculate_force(test_args.relative_permeability, test_args.first_wire_current,
            test_args.second_wire_current, test_args.length, distance)
    with raises(TypeError):
        force_law.calculate_force(test_args.relative_permeability, test_args.first_wire_current,
            test_args.second_wire_current, test_args.length, 100)
