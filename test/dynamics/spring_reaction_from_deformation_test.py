from collections import namedtuple
from pytest import approx, fixture, raises
from sympy.vector import CoordSys3D
from symplyphysics import (
    units, convert_to, SI, errors
)
from symplyphysics.core.vectors.vectors import Vector
from symplyphysics.laws.dynamics import spring_reaction_from_deformation as spring_law


@fixture
def test_args():
    C = CoordSys3D("C")
    k = units.Quantity('k')
    SI.set_quantity_dimension(k, units.force / units.length)
    SI.set_quantity_scale_factor(k, 0.1 * units.newton / units.meter)
    df_x = units.Quantity('df_x')
    SI.set_quantity_dimension(df_x, units.length)
    SI.set_quantity_scale_factor(df_x, 3 * units.meter)
    df_y = units.Quantity('df_y')
    SI.set_quantity_dimension(df_y, units.length)
    SI.set_quantity_scale_factor(df_y, 1 * units.meter)
    df = Vector([df_x, df_y], C)

    Args = namedtuple('Args', ['C', 'k', 'df'])
    return Args(C=C, k=k, df=df)

def test_basic_force(test_args):
    result = spring_law.calculate_force(test_args.k, test_args.df)
    assert SI.get_dimension_system().equivalent_dims(result.components[0].dimension, units.force)
    assert SI.get_dimension_system().equivalent_dims(result.components[1].dimension, units.force)

    result_force_x = convert_to(result.components[0], units.newton).subs(units.newton, 1).evalf(2)
    assert result_force_x == approx(-0.3, 0.01)
    result_force_y = convert_to(result.components[1], units.newton).subs(units.newton, 1).evalf(2)
    assert result_force_y == approx(-0.1, 0.01)

def test_bad_elastic_coefficient(test_args):
    eb = units.Quantity('eb')
    SI.set_quantity_dimension(eb, units.length)
    SI.set_quantity_scale_factor(eb, 1 * units.meter)

    with raises(errors.UnitsError):
        result = spring_law.calculate_force(eb, test_args.df)

    with raises(TypeError):
        result = spring_law.calculate_force(100, test_args.df)


