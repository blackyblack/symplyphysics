from collections import namedtuple
from pytest import approx, fixture, raises

from symplyphysics import (
    units, SI, errors, Probability
)
from symplyphysics.laws.nuclear import thermal_non_leakage_probability_from_diffusion_length as non_leakage_factor

@fixture
def test_args():
    thermal_diffusion_area = units.Quantity('thermal_diffusion_area')
    SI.set_quantity_dimension(thermal_diffusion_area, units.length**2)
    # water diffusion area is 8.1 cm^2
    SI.set_quantity_scale_factor(thermal_diffusion_area, 8.1 * units.centimeter**2)
    geometric_buckling = units.Quantity('geometric_buckling')
    SI.set_quantity_dimension(geometric_buckling, 1 / units.length**2)
    # sphere with radius = 1 meter
    SI.set_quantity_scale_factor(geometric_buckling, 9.869 / units.meter**2)

    Args = namedtuple('Args', ['Lth', 'Bg'])
    return Args(Lth = thermal_diffusion_area, Bg = geometric_buckling)

def test_basic_non_leakage_factor(test_args):
    result = non_leakage_factor.calculate_probability(test_args.Lth, test_args.Bg)
    assert isinstance(result, Probability)

    assert result.value == approx(0.9921, 0.01)

def test_bad_diffusion_area(test_args):
    Lthb = units.Quantity('Lthb')
    SI.set_quantity_dimension(Lthb, units.time)
    SI.set_quantity_scale_factor(Lthb, 3 * units.second)

    with raises(errors.UnitsError):
        non_leakage_factor.calculate_probability(Lthb, test_args.Bg)

    with raises(TypeError):
        non_leakage_factor.calculate_probability(100, test_args.Bg)

def test_bad_buckling(test_args):
    Bgb = units.Quantity('Bgb')
    SI.set_quantity_dimension(Bgb, units.time)
    SI.set_quantity_scale_factor(Bgb, 3 * units.second)

    with raises(errors.UnitsError):
        non_leakage_factor.calculate_probability(test_args.Lth, Bgb)

    with raises(TypeError):
        non_leakage_factor.calculate_probability(test_args.Lth, 100)

