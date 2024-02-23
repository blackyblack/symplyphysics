from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.kinematic.damped_oscillations import (
    damping_ratio_from_decay_constant_and_undamped_frequency as decay_law,
)

# Description
## A damped oscillator is described as having the decay constant zeta = 3 1/s
## and an undamped angular frequency of 6 Hz. Then its damping ratio is equal to 0.5.

Args = namedtuple("Args", "lambda_ omega")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    lambda_ = Quantity(3.0 / units.second)
    omega = Quantity(6.0 * units.hertz)
    return Args(lambda_=lambda_, omega=omega)


def test_law(test_args: Args) -> None:
    result = decay_law.calculate_damping_ratio(test_args.lambda_, test_args.omega)
    assert_equal(result, 0.5)


def test_bad_damping_constant(test_args: Args) -> None:
    lambda_bad = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        decay_law.calculate_damping_ratio(lambda_bad, test_args.omega)
    with raises(TypeError):
        decay_law.calculate_damping_ratio(100, test_args.omega)


def test_bad_frequency(test_args: Args) -> None:
    omega_bad = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        decay_law.calculate_damping_ratio(test_args.lambda_, omega_bad)
    with raises(TypeError):
        decay_law.calculate_damping_ratio(test_args.lambda_, 100)
