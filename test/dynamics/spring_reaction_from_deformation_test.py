from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    errors,
    units,
    convert_to,
    Quantity,
    SI,
)
from symplyphysics.core.vectors.vectors import QuantityVector
from symplyphysics.laws.dynamics import spring_reaction_from_deformation as spring_law


@fixture(name="test_args")
def test_args_fixture():
    k = Quantity(0.1 * units.newton / units.meter)
    df_x = Quantity(3 * units.meter)
    df_y = Quantity(1 * units.meter)
    df = QuantityVector([df_x, df_y])
    Args = namedtuple("Args", ["k", "df"])
    return Args(k=k, df=df)


def test_basic_force(test_args):
    result = spring_law.calculate_force(test_args.k, test_args.df)
    result_quantities = result.to_quantities()
    assert SI.get_dimension_system().equivalent_dims(result_quantities[0].dimension, units.force)
    assert SI.get_dimension_system().equivalent_dims(result_quantities[1].dimension, units.force)
    result_force_x = convert_to(result_quantities[0], units.newton).evalf(2)
    assert result_force_x == approx(-0.3, 0.01)
    result_force_y = convert_to(result_quantities[1], units.newton).evalf(2)
    assert result_force_y == approx(-0.1, 0.01)


def test_bad_elastic_coefficient(test_args):
    eb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        spring_law.calculate_force(eb, test_args.df)
    with raises(TypeError):
        spring_law.calculate_force(100, test_args.df)


def test_bad_deformation(test_args):
    db = Quantity(1 * units.coulomb)
    vb = QuantityVector([db])
    with raises(errors.UnitsError):
        spring_law.calculate_force(test_args.k, vb)
    with raises(TypeError):
        spring_law.calculate_force(test_args.k, 100)
