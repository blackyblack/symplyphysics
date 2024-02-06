from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_approx,
    errors,
    units,
    convert_to,
    Quantity,
    SI,
)
from symplyphysics.laws.waves import photon_energy_is_proportional_to_frequency as planck_law

# Description
## Assert we have ultraviolet radiation with frequency of 3e16 Hz.
## With online calculator
## https://www.center-pss.ru/math/raschet-energii-fotona.htm
## we obtain energy of single photone equal to 1.9878528e-17 Joule.


@fixture(name="test_args")
def test_args_fixture():
    frequency = Quantity(3e16 * units.hertz)
    Args = namedtuple("Args", ["frequency"])
    return Args(frequency=frequency)


def test_basic_energy(test_args):
    result = planck_law.calculate_energy(test_args.frequency)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.energy)
    result_energy = convert_to(result, units.joule).evalf(6)
    assert_approx(result_energy, 1.9878528e-17)


def test_bad_frequency():
    fb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        planck_law.calculate_energy(fb)
    with raises(TypeError):
        planck_law.calculate_energy(100)
