from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, units, Quantity, errors)
from symplyphysics.laws.relativistic import relativistic_momentum as momentum_law

# Description
## Let the mass be 0.1 gram, the speed is 3e6 meter per second.
## Then momentum will be equal to 300 [kilogram * meter / second].
## https://www.wolframalpha.com/input?i=relativistic+momentum+calculator&assumption=%7B%22FS%22%7D+-%3E+%7B%7B%22RelativisticMomentum%22%2C+%22p%22%7D%2C+%7B%22RelativisticMomentum%22%2C+%22v%22%7D%2C+%7B%22RelativisticMomentum%22%2C+%22m%22%7D%7D&assumption=%7B%22F%22%2C+%22RelativisticMomentum%22%2C+%22m%22%7D+-%3E%220.1+g%22&assumption=%7B%22F%22%2C+%22RelativisticMomentum%22%2C+%22v%22%7D+-%3E%223×10%5E6+m%2Fs%22


@fixture(name="test_args")
def test_args_fixture():
    mass = Quantity(0.1 * units.gram)
    velocity = Quantity(3e6 * (units.meter / units.second))

    Args = namedtuple("Args", ["mass", "velocity"])
    return Args(mass=mass, velocity=velocity)


def test_basic_momentum(test_args):
    result = momentum_law.calculate_momentum(test_args.mass, test_args.velocity)
    assert_equal(result, 300 * units.kilogram * units.meter / units.second)


def test_bad_mass(test_args):
    mass = Quantity(1 * units.joule)
    with raises(errors.UnitsError):
        momentum_law.calculate_momentum(mass, test_args.velocity)
    with raises(TypeError):
        momentum_law.calculate_momentum(100, test_args.velocity)


def test_bad_velocity(test_args):
    velocity = Quantity(1 * units.joule)
    with raises(errors.UnitsError):
        momentum_law.calculate_momentum(test_args.mass, velocity)
    with raises(TypeError):
        momentum_law.calculate_momentum(test_args.mass, 100)
