from collections import namedtuple
from pytest import approx, fixture, raises

from symplyphysics import (
    units, convert_to, SI, errors, Probability
)
from symplyphysics.laws.nuclear import thermal_utilisation_factor_from_macroscopic_absorption_cross_sections as utilisation_factor

@fixture
def test_args():
    macro_abs_fuel_cross_section = units.Quantity('macro_abs_fuel_cross_section')
    SI.set_quantity_dimension(macro_abs_fuel_cross_section, 1 / units.length)
    SI.set_quantity_scale_factor(macro_abs_fuel_cross_section, 0.2028 / units.centimeter)
    macro_abs_total_cross_section = units.Quantity('macro_abs_total_cross_section')
    SI.set_quantity_dimension(macro_abs_total_cross_section, 1 / units.length)
    SI.set_quantity_scale_factor(macro_abs_total_cross_section, 0.2356 / units.centimeter)

    Args = namedtuple('Args', ['Saf', 'Sat'])
    return Args(Saf = macro_abs_fuel_cross_section, Sat = macro_abs_total_cross_section)

def test_basic_utilisation_factor(test_args):
    result = utilisation_factor.calculate_utilisation_factor(test_args.Saf, test_args.Sat)
    assert isinstance(result, Probability)

    assert result.value == approx(0.861, 0.01)

def test_bad_macroscopic_cross_section(test_args):
    Safb = units.Quantity('Safb')
    SI.set_quantity_dimension(Safb, units.length)
    SI.set_quantity_scale_factor(Safb, 3 * units.meter)

    with raises(errors.UnitsError):
        utilisation_factor.calculate_utilisation_factor(Safb, test_args.Sat)

    with raises(TypeError):
        utilisation_factor.calculate_utilisation_factor(100, test_args.Sat)

    with raises(errors.UnitsError):
        utilisation_factor.calculate_utilisation_factor(test_args.Saf, Safb)

    with raises(TypeError):
        utilisation_factor.calculate_utilisation_factor(test_args.Saf, 100)

    Saftb = units.Quantity('Saftb')
    SI.set_quantity_dimension(Saftb, 1 / units.length)
    SI.set_quantity_scale_factor(Saftb, test_args.Sat.scale_factor + 1)
    with raises(ValueError):
        utilisation_factor.calculate_utilisation_factor(Saftb, test_args.Sat)
