from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_approx,
    errors,
    units,
    convert_to,
    Quantity,
    SI,
)
from symplyphysics.core.vectors.vectors import QuantityVector
from symplyphysics.laws.dynamics.vector import spring_reaction_from_deformation as spring_law


@fixture(name="test_args")
def test_args_fixture():
    k = Quantity(0.1 * units.newton / units.meter)
    d_x = Quantity(3 * units.meter)
    d_y = Quantity(1 * units.meter)
    d = QuantityVector([d_x, d_y])
    f_x = Quantity(-0.3 * units.newton)
    f_y = Quantity(-0.1 * units.newton)
    f = QuantityVector([f_x, f_y])
    Args = namedtuple("Args", ["k", "d", "f"])
    return Args(k=k, d=d, f=f)


def test_basic_force(test_args):
    result = spring_law.calculate_force(test_args.k, test_args.d)
    assert SI.get_dimension_system().equivalent_dims(result.components[0].dimension, units.force)
    assert SI.get_dimension_system().equivalent_dims(result.components[1].dimension, units.force)
    result_force_x = convert_to(result.components[0], units.newton).evalf(2)
    assert_approx(result_force_x, -0.3)
    result_force_y = convert_to(result.components[1], units.newton).evalf(2)
    assert_approx(result_force_y, -0.1)


def test_basic_deformation(test_args):
    result = spring_law.calculate_deformation(test_args.k, test_args.f)
    assert SI.get_dimension_system().equivalent_dims(result.components[0].dimension, units.length)
    assert SI.get_dimension_system().equivalent_dims(result.components[1].dimension, units.length)
    result_deformation_x = convert_to(result.components[0], units.meter).evalf(2)
    assert_approx(result_deformation_x, 3)
    result_deformation_y = convert_to(result.components[1], units.meter).evalf(2)
    assert_approx(result_deformation_y, 1)


def test_bad_elastic_coefficient(test_args):
    eb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        spring_law.calculate_force(eb, test_args.d)
    with raises(TypeError):
        spring_law.calculate_force(100, test_args.d)
    with raises(errors.UnitsError):
        spring_law.calculate_deformation(eb, test_args.f)
    with raises(TypeError):
        spring_law.calculate_deformation(100, test_args.f)


def test_bad_deformation(test_args):
    db = Quantity(1 * units.coulomb)
    vb = QuantityVector([db])
    with raises(errors.UnitsError):
        spring_law.calculate_force(test_args.k, vb)
    with raises(TypeError):
        spring_law.calculate_force(test_args.k, 100)


def test_bad_force(test_args):
    db = Quantity(1 * units.coulomb)
    vb = QuantityVector([db])
    with raises(errors.UnitsError):
        spring_law.calculate_deformation(test_args.k, vb)
    with raises(TypeError):
        spring_law.calculate_deformation(test_args.k, 100)
