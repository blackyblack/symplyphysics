from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.dynamics import (
    period_of_torsion_pendulum_from_rotational_inertia as pendulum_law,)

# Description
## The period of a torsion pendulum with rotational inertia I = 2.0e-4 kg*m**2 and torsion constant
## kappa = 7.1e-2 is 0.3335 s.


@fixture(name="test_args")
def test_args_fixture():
    I = Quantity(2.0e-4 * units.kilogram * units.meter**2)
    kappa = Quantity(7.1e-2 * units.newton * units.meter)
    Args = namedtuple("Args", "I kappa")
    return Args(I=I, kappa=kappa)


def test_law(test_args):
    result = pendulum_law.calculate_period(test_args.I, test_args.kappa)
    assert_equal(result, 0.3335 * units.second)


def test_bad_rotational_inertia(test_args):
    Ib = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        pendulum_law.calculate_period(Ib, test_args.kappa)
    with raises(TypeError):
        pendulum_law.calculate_period(100, test_args.kappa)


def test_bad_torsion_constant(test_args):
    kappab = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        pendulum_law.calculate_period(test_args.I, kappab)
    with raises(TypeError):
        pendulum_law.calculate_period(test_args.I, 100)
