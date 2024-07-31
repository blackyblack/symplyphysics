from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.definitions import boltzmann_factor_via_state_energy_and_temperature

# Description
## For a state with energy E = 1e-20 J of an ensemble with equilibrium temperature T = 200 K
## the Boltzmann factor is 0.0267.

Args = namedtuple("Args", "e t")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    e = Quantity(1e-20 * units.joule)
    t = Quantity(200 * units.kelvin)
    return Args(e=e, t=t)


def test_law(test_args: Args) -> None:
    result = boltzmann_factor_via_state_energy_and_temperature.calculate_boltzmann_factor(test_args.e, test_args.t)
    assert_equal(result, 0.0267, tolerance=2e-3)


def test_bad_energy(test_args: Args) -> None:
    eb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        boltzmann_factor_via_state_energy_and_temperature.calculate_boltzmann_factor(eb, test_args.t)
    with raises(TypeError):
        boltzmann_factor_via_state_energy_and_temperature.calculate_boltzmann_factor(100, test_args.t)


def test_bad_temperature(test_args: Args) -> None:
    tb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        boltzmann_factor_via_state_energy_and_temperature.calculate_boltzmann_factor(test_args.e, tb)
    with raises(TypeError):
        boltzmann_factor_via_state_energy_and_temperature.calculate_boltzmann_factor(test_args.e, 100)
