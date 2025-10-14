from collections import namedtuple
from pytest import fixture, raises
from sympy import cos, pi, sin
from symplyphysics import units, Quantity, errors
from symplyphysics.definitions.vector import net_force_vector_is_sum_of_forces as forces_law

from symplyphysics.core.experimental.coordinate_systems import CARTESIAN, QuantityCoordinateVector
from symplyphysics.core.experimental.approx import assert_equal_vectors

Args = namedtuple("Args", ["F1", "F2"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    force1 = Quantity(10 * units.newton)
    force2 = Quantity(20 * units.newton)
    F1 = QuantityCoordinateVector([
        force1 * cos(pi / 3),
        force1 * sin(pi / 3),
        0,
    ], CARTESIAN)
    F2 = QuantityCoordinateVector([
        force2,
        0,
        0,
    ], CARTESIAN)
    return Args(F1=F1, F2=F2)


def test_basic_superposition(test_args: Args) -> None:
    result = forces_law.calculate_resultant_force([test_args.F1, test_args.F2])

    expected = QuantityCoordinateVector([
        25 * units.newton,
        8.66 * units.newton,
        0,
    ], CARTESIAN)

    assert_equal_vectors(result, expected)


def test_three_forces_array(test_args: Args) -> None:
    F3 = QuantityCoordinateVector([
        Quantity(0, dimension=units.force),
        Quantity(-5 * units.newton),
        0,
    ], CARTESIAN)

    result = forces_law.calculate_resultant_force([test_args.F1, test_args.F2, F3])

    expected = QuantityCoordinateVector([
        25 * units.newton,
        3.66 * units.newton,
        0,
    ], CARTESIAN)

    assert_equal_vectors(result, expected)


def test_bad_force(test_args: Args) -> None:
    Fb = QuantityCoordinateVector([
        Quantity(1 * units.meter),
        0,
        0,
    ], CARTESIAN)

    with raises(errors.UnitsError):
        forces_law.calculate_resultant_force([Fb, test_args.F2])
    with raises(TypeError):
        forces_law.calculate_resultant_force([100, test_args.F2])
    with raises(errors.UnitsError):
        forces_law.calculate_resultant_force([test_args.F1, Fb])
    with raises(TypeError):
        forces_law.calculate_resultant_force([test_args.F1, 100])
    with raises(errors.UnitsError):
        forces_law.calculate_resultant_force([Fb, Fb])
    with raises(TypeError):
        forces_law.calculate_resultant_force([100, 100])
    with raises(TypeError):
        forces_law.calculate_resultant_force(test_args.F1)

    assert forces_law.calculate_resultant_force([]) == 0
