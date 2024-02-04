from collections import namedtuple
from pytest import fixture, raises
from sympy import cos, pi, sin
from symplyphysics import (
    assert_approx,
    units,
    convert_to,
    Quantity,
    SI,
)
from symplyphysics.core import errors
from symplyphysics.core.vectors.vectors import QuantityVector
from symplyphysics.definitions.vector import superposition_of_forces_is_sum as forces_law


@fixture(name="test_args")
def test_args_fixture():
    force1 = Quantity(10 * units.newton)
    force2 = Quantity(20 * units.newton)
    F1 = QuantityVector([force1 * cos(pi / 3), force1 * sin(pi / 3)])
    F2 = QuantityVector([force2, 0])
    Args = namedtuple("Args", ["F1", "F2"])
    return Args(F1=F1, F2=F2)


def test_basic_superposition(test_args):
    result = forces_law.calculate_resultant_force([test_args.F1, test_args.F2])
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.force)
    result_force_x = convert_to(result.components[0], units.newton).evalf(3)
    assert_approx(result_force_x, 25)
    result_force_y = convert_to(result.components[1], units.newton).evalf(3)
    assert_approx(result_force_y, 8.66)
    assert len(result.components) == 2


def test_three_forces_array(test_args):
    F3 = QuantityVector([Quantity(0, dimension=units.force), Quantity(-5 * units.newton)])
    result = forces_law.calculate_resultant_force([test_args.F1, test_args.F2, F3])
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.force)
    result_force_x = convert_to(result.components[0], units.newton).evalf(3)
    assert_approx(result_force_x, 25)
    result_force_y = convert_to(result.components[1], units.newton).evalf(3)
    assert_approx(result_force_y, 3.66)


def test_bad_force(test_args):
    Fb = QuantityVector([Quantity(1 * units.meter)])
    with raises(errors.UnitsError):
        forces_law.calculate_resultant_force([Fb, test_args.F2])
    with raises(AttributeError):
        forces_law.calculate_resultant_force([100, test_args.F2])
    with raises(errors.UnitsError):
        forces_law.calculate_resultant_force([test_args.F1, Fb])
    with raises(AttributeError):
        forces_law.calculate_resultant_force([test_args.F1, 100])
    with raises(errors.UnitsError):
        forces_law.calculate_resultant_force([Fb, Fb])
    with raises(AttributeError):
        forces_law.calculate_resultant_force([100, 100])
    with raises(TypeError):
        forces_law.calculate_resultant_force(test_args.F1)
