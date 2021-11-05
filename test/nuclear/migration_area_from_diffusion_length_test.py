from collections import namedtuple
from pytest import approx, fixture, raises

from symplyphysics import (
    units, convert_to, SI, errors
)
from symplyphysics.laws.nuclear import migration_area_from_diffusion_length as migration_area

@fixture
def test_args():
    diffusion_area = units.Quantity('diffusion_area')
    SI.set_quantity_dimension(diffusion_area, units.length**2)
    # water diffusion area = 39.9 cm^2
    SI.set_quantity_scale_factor(diffusion_area, 39.9 * units.centimeter**2)
    neutron_fermi_age = units.Quantity('neutron_fermi_age')
    SI.set_quantity_dimension(neutron_fermi_age, units.length**2)
    # light water Fermi age = 38.8 cm^2
    SI.set_quantity_scale_factor(neutron_fermi_age, 38.8 * units.centimeter**2)

    Args = namedtuple('Args', ['Ld', 'th'])
    return Args(Ld = diffusion_area, th = neutron_fermi_age)

def test_basic_migration_area(test_args):
    result = migration_area.calculate_migration_area(test_args.Ld, test_args.th)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.length**2)

    result_area = convert_to(result, units.centimeter**2).subs(units.centimeter, 1).evalf(2)
    # water migration area = 78.8 cm^2
    assert result_area == approx(78.8, 0.01)

def test_bad_diffusion_area(test_args):
    Ldb = units.Quantity('Ldb')
    SI.set_quantity_dimension(Ldb, units.time)
    SI.set_quantity_scale_factor(Ldb, 3 * units.second)

    with raises(errors.UnitsError):
        migration_area.calculate_migration_area(Ldb, test_args.th)

    with raises(TypeError):
        migration_area.calculate_migration_area(100, test_args.th)

def test_bad_macroscopic_cross_section(test_args):
    thb = units.Quantity('thb')
    SI.set_quantity_dimension(thb, units.time)
    SI.set_quantity_scale_factor(thb, 3 * units.second)

    with raises(errors.UnitsError):
        migration_area.calculate_migration_area(test_args.Ld, thb)

    with raises(TypeError):
        migration_area.calculate_migration_area(test_args.Ld, 100)
