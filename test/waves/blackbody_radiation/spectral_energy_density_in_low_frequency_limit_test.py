from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    prefixes,
    Quantity,
)
from symplyphysics.laws.waves.blackbody_radiation import spectral_energy_density_in_low_frequency_limit as rj_law

# Description
## The spectral energy density of electromagnetic radiation emitted from the Sun (assuming it is a black
## body radiator, the surface temperature of the Sun is about T = 5800 K) for green light (nu = 550 THz) is
## u_nu = 141 keV / (m**3 * Hz), assuming Rayleigh-Jeans law approximation.

Args = namedtuple("Args", "nu t")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    nu = Quantity(550 * prefixes.tera * units.hertz)
    t = Quantity(5800 * units.kelvin)
    return Args(nu=nu, t=t)


def test_law(test_args: Args) -> None:
    result = rj_law.calculate_spectral_energy_density(test_args.nu, test_args.t)
    assert_equal(result, 141 * prefixes.kilo * units.electronvolt / (units.meter**3 * units.hertz))


def test_bad_frequency(test_args: Args) -> None:
    nub = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        rj_law.calculate_spectral_energy_density(nub, test_args.t)
    with raises(TypeError):
        rj_law.calculate_spectral_energy_density(100, test_args.t)


def test_bad_temperature(test_args: Args) -> None:
    tb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        rj_law.calculate_spectral_energy_density(test_args.nu, tb)
    with raises(TypeError):
        rj_law.calculate_spectral_energy_density(test_args.nu, 100)
