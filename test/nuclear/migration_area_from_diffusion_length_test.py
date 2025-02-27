from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.nuclear import migration_area_from_diffusion_length as migration_area

Args = namedtuple("Args", ["Ld", "th"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    # water diffusion area = 39.9 cm^2
    diffusion_area = Quantity(39.9 * units.centimeter**2)
    # light water Fermi age = 38.8 cm^2
    neutron_fermi_age = Quantity(38.8 * units.centimeter**2)
    return Args(Ld=diffusion_area, th=neutron_fermi_age)


def test_basic_migration_area(test_args: Args) -> None:
    result = migration_area.calculate_migration_area(test_args.Ld, test_args.th)
    # water migration area = 78.8 cm^2
    assert_equal(result, 78.8 * units.centimeter**2, relative_tolerance=0.01)


def test_bad_diffusion_area(test_args: Args) -> None:
    Lb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        migration_area.calculate_migration_area(Lb, test_args.th)
    with raises(TypeError):
        migration_area.calculate_migration_area(100, test_args.th)


def test_bad_fermi_age(test_args: Args) -> None:
    thb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        migration_area.calculate_migration_area(test_args.Ld, thb)
    with raises(TypeError):
        migration_area.calculate_migration_area(test_args.Ld, 100)
