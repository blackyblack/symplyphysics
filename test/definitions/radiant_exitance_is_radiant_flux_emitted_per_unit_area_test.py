from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    prefixes,
    Quantity,
)
from symplyphysics.definitions import (
    radiant_exitance_is_radiant_flux_emitted_per_unit_area as radiant_exitance_def,)

# Description
## The radiant exitance of a surface of area A = 0.1 mm**2 which emits radiant flux Phi_e = 1 mW
## amounts to 10 kW/m**2.

Args = namedtuple("Args", "flux area")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    flux = Quantity(1 * prefixes.milli * units.watt)
    area = Quantity(0.1 * units.millimeter**2)
    return Args(flux=flux, area=area)


def test_law(test_args: Args) -> None:
    result = radiant_exitance_def.calculate_radiant_exitance(test_args.flux, test_args.area)
    assert_equal(result, 10 * prefixes.kilo * units.watt / units.meter**2)


def test_bad_flux(test_args: Args) -> None:
    bad_flux = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        radiant_exitance_def.calculate_radiant_exitance(bad_flux, test_args.area)
    with raises(TypeError):
        radiant_exitance_def.calculate_radiant_exitance(100, test_args.area)


def test_bad_area(test_args: Args) -> None:
    bad_area = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        radiant_exitance_def.calculate_radiant_exitance(test_args.flux, bad_area)
    with raises(TypeError):
        radiant_exitance_def.calculate_radiant_exitance(test_args.flux, 100)
