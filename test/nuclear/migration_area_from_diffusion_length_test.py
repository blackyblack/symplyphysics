from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    errors,
    units,
    convert_to,
    Quantity,
    SI,
)
from symplyphysics.laws.nuclear import migration_area_from_diffusion_length as migration_area


@fixture
def test_args():
    # water diffusion area = 39.9 cm^2
    diffusion_area = Quantity(39.9 * units.centimeter**2)
    # light water Fermi age = 38.8 cm^2
    neutron_fermi_age = Quantity(38.8 * units.centimeter**2)
    Args = namedtuple("Args", ["Ld", "th"])
    return Args(Ld=diffusion_area, th=neutron_fermi_age)


def test_basic_migration_area(test_args):
    result = migration_area.calculate_migration_area(test_args.Ld, test_args.th)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.length**2)
    result_area = convert_to(result, units.centimeter**2).subs(units.centimeter, 1).evalf(2)
    # water migration area = 78.8 cm^2
    assert result_area == approx(78.8, 0.01)


def test_bad_diffusion_area(test_args):
    Ldb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        migration_area.calculate_migration_area(Ldb, test_args.th)
    with raises(TypeError):
        migration_area.calculate_migration_area(100, test_args.th)


def test_bad_fermi_age(test_args):
    thb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        migration_area.calculate_migration_area(test_args.Ld, thb)
    with raises(TypeError):
        migration_area.calculate_migration_area(test_args.Ld, 100)
