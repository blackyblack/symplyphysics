from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.thermodynamics.maxwell_boltzmann_statistics import energy_distribution

# Description
## The value of the energy distribution function for an ensemble of Argon particles at
## equilibrium temperature T = 300 K around energy value E = 1.24 eV is f(E) = 2.78 1/J.

Args = namedtuple("Args", "e t")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    e = Quantity(1.24 * units.electronvolt)
    t = Quantity(300 * units.kelvin)
    return Args(e=e, t=t)


def test_law(test_args: Args) -> None:
    result = energy_distribution.calculate_energy_distribution_function(test_args.e, test_args.t)
    assert_equal(result, 2.78 / units.joule, tolerance=2e-3)


def test_bad_energy(test_args: Args) -> None:
    eb = Quantity(1.0 * units.meter)
    with raises(errors.UnitsError):
        energy_distribution.calculate_energy_distribution_function(eb, test_args.t)
    with raises(TypeError):
        energy_distribution.calculate_energy_distribution_function(100, test_args.t)


def test_bad_temperature(test_args: Args) -> None:
    tb = Quantity(1.0 * units.meter)
    with raises(errors.UnitsError):
        energy_distribution.calculate_energy_distribution_function(test_args.e, tb)
    with raises(TypeError):
        energy_distribution.calculate_energy_distribution_function(test_args.e, 100)
