from collections import namedtuple
from pytest import approx, fixture, raises

from symplyphysics import (
    units, convert_to, SI, errors
)
from symplyphysics.laws.nuclear.buckling import material_buckling_from_macroscopic_fission_cross_section_diffusion_coefficient as material_buckling

@fixture
def test_args():
    neutrons_per_fission = 2.6
    macro_fission_cross_section = units.Quantity('macro_fission_cross_section')
    SI.set_quantity_dimension(macro_fission_cross_section, 1 / units.length)
    SI.set_quantity_scale_factor(macro_fission_cross_section, 1.482 / units.centimeter)
    macro_abs_cross_section = units.Quantity('macro_abs_cross_section')
    SI.set_quantity_dimension(macro_abs_cross_section, 1 / units.length)
    SI.set_quantity_scale_factor(macro_abs_cross_section, 3.108 / units.centimeter)
    diffusion_coefficient = units.Quantity('diffusion_coefficient')
    SI.set_quantity_dimension(diffusion_coefficient, units.length)
    SI.set_quantity_scale_factor(diffusion_coefficient, 31.782 * units.meter)

    Args = namedtuple('Args', ['v', 'Sf', 'Sa', 'D'])
    return Args(v = neutrons_per_fission, Sf = macro_fission_cross_section, Sa = macro_abs_cross_section,
        D = diffusion_coefficient)

def test_basic_buckling(test_args):
    result = material_buckling.calculate_buckling(test_args.v, test_args.Sf, test_args.Sa, test_args.D)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, 1 / units.length**2)

    result_buckling = convert_to(result, 1 / units.meter**2).subs(units.meter, 1).evalf(2)
    assert result_buckling == approx(2.345, 0.01)

def test_bad_macroscopic_cross_section(test_args):
    Sb = units.Quantity('Sb')
    SI.set_quantity_dimension(Sb, units.time)
    SI.set_quantity_scale_factor(Sb, 3 * units.second)

    with raises(errors.UnitsError):
        material_buckling.calculate_buckling(test_args.v, Sb, test_args.Sa, test_args.D)

    with raises(TypeError):
        material_buckling.calculate_buckling(test_args.v, 100, test_args.Sa, test_args.D)

    with raises(errors.UnitsError):
        material_buckling.calculate_buckling(test_args.v, test_args.Sf, Sb, test_args.D)

    with raises(TypeError):
        material_buckling.calculate_buckling(test_args.v, test_args.Sf, 100, test_args.D)

def test_bad_diffusion_coefficient(test_args):
    Db = units.Quantity('Db')
    SI.set_quantity_dimension(Db, units.time)
    SI.set_quantity_scale_factor(Db, 3 * units.second)

    with raises(errors.UnitsError):
        material_buckling.calculate_buckling(test_args.v, test_args.Sf, test_args.Sa, Db)

    with raises(TypeError):
        material_buckling.calculate_buckling(test_args.v, test_args.Sf, test_args.Sa, 100)
