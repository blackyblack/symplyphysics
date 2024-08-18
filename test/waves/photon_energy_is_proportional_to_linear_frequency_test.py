from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.waves import photon_energy_is_proportional_to_linear_frequency as planck_law

# Description
## Assert we have ultraviolet radiation with frequency of 3e16 Hz.
## With online calculator
## https://www.center-pss.ru/math/raschet-energii-fotona.htm
## we obtain energy of single photone equal to 1.9878528e-17 Joule.

Args = namedtuple("Args", ["frequency"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    frequency = Quantity(3e16 * units.hertz)
    return Args(frequency=frequency)


def test_basic_energy(test_args: Args) -> None:
    result = planck_law.calculate_energy(test_args.frequency)
    assert_equal(result, 1.9878528e-17 * units.joule)


def test_bad_frequency() -> None:
    fb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        planck_law.calculate_energy(fb)
    with raises(TypeError):
        planck_law.calculate_energy(100)
