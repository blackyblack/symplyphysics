from collections import namedtuple
from pytest import fixture, raises
from sympy import cos, pi, sin
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.core.coordinate_systems.coordinate_systems import CoordinateSystem
from symplyphysics.core.vectors.vectors import QuantityVector
from symplyphysics.laws.dynamics.vector import mechanical_work_from_force_and_move as work_law

# Description
## Force of 100N is applied to heavy object lying oh horizontal table. Force is directed 60 degrees up.
## Object slides on the table with no friction and is been moved to 3m far.
## As cosine of 60 degrees is 1/2, work should be 100 * 3 / 2 = 150 Joules.

Args = namedtuple("Args", ["F", "S"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    force_ = Quantity(100 * units.newton)
    distance_ = Quantity(3 * units.meter)
    coordinates = CoordinateSystem(CoordinateSystem.System.CARTESIAN)
    force_vector = QuantityVector([force_ * cos(pi / 3), force_ * sin(pi / 3)], coordinates)
    distance_vector = QuantityVector([distance_, 0], coordinates)
    return Args(F=force_vector, S=distance_vector)


def test_basic_work(test_args: Args) -> None:
    result = work_law.calculate_work(test_args.F, test_args.S)
    assert_equal(result, 150 * units.joule)


def test_bad_force(test_args: Args) -> None:
    Fb = Quantity(1 * units.coulomb)
    force_vector = QuantityVector([Fb])
    with raises(errors.UnitsError):
        work_law.calculate_work(force_vector, test_args.S)
    with raises(TypeError):
        work_law.calculate_work(100, test_args.S)


def test_bad_move(test_args: Args) -> None:
    Sb = Quantity(1 * units.coulomb)
    move_vector = QuantityVector([Sb])
    with raises(errors.UnitsError):
        work_law.calculate_work(test_args.F, move_vector)
    with raises(TypeError):
        work_law.calculate_work(test_args.F, 100)
