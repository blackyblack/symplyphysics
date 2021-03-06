from collections import namedtuple
from pytest import approx, fixture, raises

from symplyphysics import (
    units, convert_to, SI, errors
)
from symplyphysics.laws.nuclear import diffusion_area_from_diffusion_coefficient as diffusion_area

@fixture
def test_args():
    macro_abs_cross_section = units.Quantity('macro_abs_cross_section')
    SI.set_quantity_dimension(macro_abs_cross_section, 1 / units.length)
    # water macroscopic absorption cross-section is 0.022 cm^-1
    SI.set_quantity_scale_factor(macro_abs_cross_section, 0.022 / units.centimeter)
    diffusion_coefficient = units.Quantity('diffusion_coefficient')
    SI.set_quantity_dimension(diffusion_coefficient, units.length)
    # water diffusion coefficient is 0.142 cm
    SI.set_quantity_scale_factor(diffusion_coefficient, 0.142 * units.centimeter)

    Args = namedtuple('Args', ['Sa', 'D'])
    return Args(Sa = macro_abs_cross_section, D = diffusion_coefficient)

def test_basic_diffusion_length(test_args):
    result = diffusion_area.calculate_diffusion_area(test_args.D, test_args.Sa)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.length**2)

    result_diffusion = convert_to(result, units.centimeter**2).subs(units.centimeter, 1).evalf(2)
    assert result_diffusion == approx(2.54**2, 0.01)

def test_bad_diffusion_coefficient(test_args):
    Db = units.Quantity('Db')
    SI.set_quantity_dimension(Db, units.time)
    SI.set_quantity_scale_factor(Db, 3 * units.second)

    with raises(errors.UnitsError):
        diffusion_area.calculate_diffusion_area(Db, test_args.Sa)

    with raises(TypeError):
        diffusion_area.calculate_diffusion_area(100, test_args.Sa)

def test_bad_macroscopic_cross_section(test_args):
    Sb = units.Quantity('Sb')
    SI.set_quantity_dimension(Sb, units.time)
    SI.set_quantity_scale_factor(Sb, 3 * units.second)

    with raises(errors.UnitsError):
        diffusion_area.calculate_diffusion_area(test_args.D, Sb)

    with raises(TypeError):
        diffusion_area.calculate_diffusion_area(test_args.D, 100)
