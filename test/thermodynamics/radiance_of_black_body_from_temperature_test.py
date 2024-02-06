from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.thermodynamics import radiance_of_black_body_from_temperature as stefan_boltzmann_law

# Description. With help of online calculator at https://www.wolframalpha.com/input/?i=stefan+boltzmann+law:
## For black body with temperature 20 K radiant heat energy should be 0.009073 [W/m^2]


@fixture(name="test_args")
def test_args_fixture():
    t = Quantity(20 * units.kelvin)
    Args = namedtuple("Args", ["t"])
    return Args(t=t)


def test_basic_radiance_heat_energy(test_args):
    result = stefan_boltzmann_law.calculate_radiance(test_args.t)
    assert_equal(result, 0.009073 * units.watt / units.meter**2)


def test_bad_temperature():
    tb = Quantity(1 * units.second)
    with raises(errors.UnitsError):
        stefan_boltzmann_law.calculate_radiance(tb)
    with raises(TypeError):
        stefan_boltzmann_law.calculate_radiance(100)
