from collections import namedtuple
from pytest import approx, fixture, raises

from symplyphysics import (
    units, convert_to, SI, errors
)
from symplyphysics.laws.nuclear import neutron_diffusion_coefficient_from_scattering_cross_section as diffusion_coeff

@fixture
def test_args():
    macro_transport_cross_section = units.Quantity('macro_transport_cross_section')
    SI.set_quantity_dimension(macro_transport_cross_section, 1 / units.length)
    # carbon macroscopic transport cross-section
    SI.set_quantity_scale_factor(macro_transport_cross_section, 0.4987 / units.centimeter)

    Args = namedtuple('Args', ['S'])
    return Args(S = macro_transport_cross_section)

def test_basic_coefficient(test_args):
    result = diffusion_coeff.calculate_diffusion_coefficient(test_args.S)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.length)

    result_coefficient = convert_to(result, units.centimeter).subs(units.centimeter, 1).evalf(2)
    # carbon diffusion coefficient is 0.668 cm
    assert result_coefficient == approx(0.668, 0.01)

def test_bad_macroscopic_cross_section():
    Sb = units.Quantity('Sb')
    SI.set_quantity_dimension(Sb, units.time)
    SI.set_quantity_scale_factor(Sb, 3 * units.second)

    with raises(errors.UnitsError):
        diffusion_coeff.calculate_diffusion_coefficient(Sb)

    with raises(TypeError):
        diffusion_coeff.calculate_diffusion_coefficient(100)
