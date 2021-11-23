from collections import namedtuple
from pytest import approx, fixture, raises

from sympy.core.singleton import S
from sympy.physics.units import Dimension
from symplyphysics import (
    units, convert_to, SI, errors
)
from symplyphysics.laws.nuclear.buckling import geometric_buckling_from_infinite_multiplication_factor_diffusion_area as buckling

@fixture
def test_args():
    infinite_multiplication_factor = 1.0247
    # critical reactor
    effective_multiplication_factor = 1
    diffusion_area = units.Quantity('diffusion_area')
    SI.set_quantity_dimension(diffusion_area, units.length**2)
    SI.set_quantity_scale_factor(diffusion_area, 60.117 * units.centimeter**2)

    Args = namedtuple('Args', ['k_inf', 'k_eff', 'L2'])
    return Args(k_inf = infinite_multiplication_factor, k_eff = effective_multiplication_factor, L2 = diffusion_area)

def test_basic_buckling(test_args):
    result = buckling.calculate_geometric_buckling_squared(test_args.k_inf, test_args.k_eff, test_args.L2)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.length**-2)

    result_buckling = convert_to(result, units.centimeter**-2).subs(units.centimeter, 1).evalf(4)
    assert result_buckling == approx(0.000412, 0.01)

def test_zero_buckling(test_args):
    result = buckling.calculate_geometric_buckling_squared(1, 1, test_args.L2)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, Dimension(1))

    result_buckling = convert_to(result, S.One).evalf(4)
    assert result_buckling == approx(0, 0.01)

def test_bad_diffusion_area(test_args):
    L2b = units.Quantity('L2b')
    SI.set_quantity_dimension(L2b, units.time)
    SI.set_quantity_scale_factor(L2b, 3 * units.second)

    with raises(errors.UnitsError):
        buckling.calculate_geometric_buckling_squared(test_args.k_inf, test_args.k_eff, L2b)

    with raises(TypeError):
        buckling.calculate_geometric_buckling_squared(test_args.k_inf, test_args.k_eff, 100)
