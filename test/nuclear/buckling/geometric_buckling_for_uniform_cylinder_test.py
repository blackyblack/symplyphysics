from collections import namedtuple
from pytest import approx, fixture, raises

from symplyphysics import (
    units, convert_to, SI, errors
)
from symplyphysics.laws.nuclear.buckling import geometric_buckling_for_uniform_cylinder as buckling

@fixture
def test_args():
    cylinder_radius = units.Quantity('cylinder_radius')
    SI.set_quantity_dimension(cylinder_radius, units.length)
    SI.set_quantity_scale_factor(cylinder_radius, 1.8 * units.meter)
    cylinder_height = units.Quantity('cylinder_height')
    SI.set_quantity_dimension(cylinder_height, units.length)
    SI.set_quantity_scale_factor(cylinder_height, 4.2 * units.meter)

    Args = namedtuple('Args', ['R', 'H'])
    return Args(R = cylinder_radius, H = cylinder_height)

def test_basic_geometric_buckling(test_args):
    result = buckling.calculate_geometric_buckling_squared(test_args.R, test_args.H)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.length**-2)

    result_geometric_buckling = convert_to(result, units.meter**-2).subs(units.meter, 1).evalf(2)
    assert result_geometric_buckling == approx(2.3447, 0.01)

def test_bad_cylinder_radius(test_args):
    Rb = units.Quantity('Rb')
    SI.set_quantity_dimension(Rb, units.time)
    SI.set_quantity_scale_factor(Rb, 3 * units.second)

    with raises(errors.UnitsError):
        buckling.calculate_geometric_buckling_squared(Rb, test_args.H)

    with raises(TypeError):
        buckling.calculate_geometric_buckling_squared(100, test_args.H)

def test_bad_cylinder_height(test_args):
    Hb = units.Quantity('Hb')
    SI.set_quantity_dimension(Hb, units.time)
    SI.set_quantity_scale_factor(Hb, 3 * units.second)

    with raises(errors.UnitsError):
        buckling.calculate_geometric_buckling_squared(test_args.R, Hb)

    with raises(TypeError):
        buckling.calculate_geometric_buckling_squared(test_args.R, 100)
