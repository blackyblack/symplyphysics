from collections import namedtuple
from pytest import fixture, raises
from sympy import cos, pi, sin
from symplyphysics import (
    assert_equal_vectors,
    units,
    Quantity,
)
from symplyphysics.core import errors
from symplyphysics.core.vectors.vectors import QuantityVector
from symplyphysics.definitions.vector import superposition_of_forces_is_sum as forces_law

Args = namedtuple("Args", ["F1", "F2"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    force1 = Quantity(10 * units.newton)
    force2 = Quantity(20 * units.newton)
    F1 = QuantityVector([force1 * cos(pi / 3), force1 * sin(pi / 3)])
    F2 = QuantityVector([force2, 0])
    return Args(F1=F1, F2=F2)


def test_basic_superposition(test_args: Args) -> None:
    result = forces_law.calculate_resultant_force([test_args.F1, test_args.F2])
    assert_equal_vectors(
        result,
        QuantityVector([25 * units.newton, 8.66 * units.newton, 0]),
    )


def test_three_forces_array(test_args: Args) -> None:
    F3 = QuantityVector([Quantity(0, dimension=units.force), Quantity(-5 * units.newton)])
    result = forces_law.calculate_resultant_force([test_args.F1, test_args.F2, F3])
    assert_equal_vectors(
        result,
        QuantityVector([25 * units.newton, 3.66 * units.newton, 0]),
    )


def test_bad_force(test_args: Args) -> None:
    Fb = QuantityVector([Quantity(1 * units.meter)])
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
    with raises(ValueError):
        forces_law.calculate_resultant_force([])
