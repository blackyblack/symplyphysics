from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, units, Quantity, errors)
from symplyphysics.laws.electricity import electric_field_of_infinite_charged_plane as intensity_law

# Description
## It is known that with a surface charge density of 2.5 [coulomb / meter^2], the electric field strength is 1.4e11 [volt / meter].
## https://www.calculatoratoz.com/en/electric-field-due-to-infinite-sheet-calculator/Calc-672

Args = namedtuple("Args", ["surface_charge_density"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    surface_charge_density = Quantity(2.5 * (units.coulomb / units.meter**2))
    return Args(surface_charge_density=surface_charge_density)


def test_basic_electric_intensity(test_args: Args) -> None:
    result = intensity_law.calculate_electric_intensity(test_args.surface_charge_density)
    assert_equal(result, 1.412e11 * units.volt / units.meter)


def test_bad_surface_charge_density() -> None:
    surface_charge_density = Quantity(1 * units.joule)
    with raises(errors.UnitsError):
        intensity_law.calculate_electric_intensity(surface_charge_density)
    with raises(TypeError):
        intensity_law.calculate_electric_intensity(100)
