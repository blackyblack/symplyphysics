from collections import namedtuple
from pytest import approx, fixture, raises

from symplyphysics import (
    units, convert_to, SI, errors
)
from symplyphysics.laws.nuclear import macroscopic_cross_section_from_microscopic_cross_section as macro_cs

@fixture
def test_args():
    microscopic_cross_section = units.Quantity('microscopic_cross_section')
    SI.set_quantity_dimension(microscopic_cross_section, units.length**2)
    # carbon microscopic cross-section is 5.01 barn
    SI.set_quantity_scale_factor(microscopic_cross_section, 5.01E-24 * units.centimeter**2)
    atomic_number_density = units.Quantity('atomic_number_density')
    SI.set_quantity_dimension(atomic_number_density, 1 / units.length**3)
    # NC = 1 * 2.75*10^22 atoms of carbon/cm^3
    SI.set_quantity_scale_factor(atomic_number_density, 2.75E+22 / units.centimeter**3)

    Args = namedtuple('Args', ['b', 'N'])
    return Args(b = microscopic_cross_section, N = atomic_number_density)

def test_basic_cross_section(test_args):
    result = macro_cs.calculate_cross_section(test_args.b, test_args.N)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.length**-1)

    result_cross_section = convert_to(result, units.centimeter**-1).subs(units.centimeter, 1).evalf(2)
    assert result_cross_section == approx(0.14, 0.1)

def test_bad_microscopic_cross_section(test_args):
    bb = units.Quantity('bb')
    SI.set_quantity_dimension(bb, units.length)
    SI.set_quantity_scale_factor(bb, 3 * units.meter)

    with raises(errors.UnitsError):
        macro_cs.calculate_cross_section(bb, test_args.N)

    with raises(TypeError):
        macro_cs.calculate_cross_section(100, test_args.N)

def test_bad_atomic_number_density(test_args):
    Nb = units.Quantity('Nb')
    SI.set_quantity_dimension(Nb, units.length)
    SI.set_quantity_scale_factor(Nb, 3 * units.meter)

    with raises(errors.UnitsError):
        macro_cs.calculate_cross_section(test_args.b, Nb)

    with raises(TypeError):
        macro_cs.calculate_cross_section(test_args.b, 100)
