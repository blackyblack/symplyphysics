from collections import namedtuple
from pytest import approx, fixture, raises

from symplyphysics import (
    units, SI, errors
)
from symplyphysics.laws.nuclear import resonance_escape_probability_from_resonance_absorption_integral as resonance_escape

@fixture
def test_args():
    atomic_number_density_abs = units.Quantity('atomic_number_density_abs')
    SI.set_quantity_dimension(atomic_number_density_abs, 1 / units.length ** 3)
    # U-238 atomic number density = 0.984 atoms per (barn * cm)
    SI.set_quantity_scale_factor(atomic_number_density_abs, 0.984E+24 / units.centimeter ** 3)
    resonance_integral = units.Quantity('resonance_integral')
    SI.set_quantity_dimension(resonance_integral, units.length ** 2)
    # UO2 fuel rod with diameter = 1 cm has effective resonance integral = 20.5 barns
    SI.set_quantity_scale_factor(resonance_integral, 20.5E-24 * units.centimeter ** 2)
    # U-238 with carbon moderator average lethargy decrement = 0.1573
    average_lethargy_change = 0.1573
    macro_scatter_cross_section = units.Quantity('macro_scatter_cross_section')
    SI.set_quantity_dimension(macro_scatter_cross_section, 1 / units.length)
    # Carbon macroscopic cross-section = 2820 cm^-1
    SI.set_quantity_scale_factor(macro_scatter_cross_section, 2820 / units.centimeter)

    Args = namedtuple('Args', ['Na', 'Ieff', 'Let', 'Ss'])
    return Args(Na = atomic_number_density_abs, Ieff = resonance_integral, Let = average_lethargy_change, Ss = macro_scatter_cross_section)

def test_basic_resonance_escape_factor(test_args):
    result = resonance_escape.calculate_resonance_escape_probability(test_args.Na, test_args.Ieff, test_args.Let, test_args.Ss)
    assert result.value == approx(0.955, 0.01)

def test_bad_atomic_number_density(test_args):
    Nab = units.Quantity('Nab')
    SI.set_quantity_dimension(Nab, units.time)
    SI.set_quantity_scale_factor(Nab, 3 * units.second)

    with raises(errors.UnitsError):
        resonance_escape.calculate_resonance_escape_probability(Nab, test_args.Ieff, test_args.Let, test_args.Ss)

    with raises(TypeError):
        resonance_escape.calculate_resonance_escape_probability(100, test_args.Ieff, test_args.Let, test_args.Ss)
