from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.dynamics import friction_force_from_normal_force as friction_force_law

#Description. According to online friction calculator (https://www.omnicalculator.com/physics/friction) normal reaction of 5 Newtons with friction factor of 0.001 will cause 0.005 newtons of friction force.


@fixture(name="test_args")
def test_args_fixture():
    mu = 0.001
    N = Quantity(5 * units.newton)
    Args = namedtuple("Args", ["mu", "N"])
    return Args(mu=mu, N=N)


def test_basic_friction_force(test_args):
    result = friction_force_law.calculate_friction_force(test_args.mu, test_args.N)
    assert_equal(result, 0.005 * units.newton)


def test_bad_reaction(test_args):
    Nb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        friction_force_law.calculate_friction_force(test_args.mu, Nb)
    with raises(TypeError):
        friction_force_law.calculate_friction_force(test_args.mu, 100)
