from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, units, Quantity, errors)
from symplyphysics.laws.condensed_matter import current_density_in_thermionic_emission_per_richardson as emission_law

# Description
## It is known that for strontium, the current density of thermionic emission at
## a temperature of 300 kelvin is equal to 0.33e-32 [A / m^2].
## https://saecanet.com/Calculation_data/000034_000170_richardson-dushman.html

Args = namedtuple("Args", ["thermodynamic_work", "temperature"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    thermodynamic_work = Quantity(2.59 * units.electronvolt)
    temperature = Quantity(300 * units.kelvin)

    return Args(thermodynamic_work=thermodynamic_work, temperature=temperature)


def test_basic_thermionic_current(test_args: Args) -> None:
    result = emission_law.calculate_current(test_args.thermodynamic_work, test_args.temperature)
    assert_equal(result, 3.342e-33 * units.ampere / units.meter**2)


def test_bad_thermodynamic_work(test_args: Args) -> None:
    thermodynamic_work = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        emission_law.calculate_current(thermodynamic_work, test_args.temperature)
    with raises(TypeError):
        emission_law.calculate_current(100, test_args.temperature)


def test_bad_temperature(test_args: Args) -> None:
    temperature = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        emission_law.calculate_current(test_args.thermodynamic_work, temperature)
    with raises(TypeError):
        emission_law.calculate_current(test_args.thermodynamic_work, 100)
