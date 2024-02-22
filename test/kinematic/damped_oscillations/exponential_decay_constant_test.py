from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.kinematic.damped_oscillations import (
    exponential_decay_constant as decay_law,
)

# Description
## A damped oscillator has a damping ratio of 0.5 and its undamped angular frequency
## is 6 Hz. Then the exponential decay constant of the oscillator is 3 1/s.

Args = namedtuple("Args", "zeta omega")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    zeta = 0.5
    omega = Quantity(6.0 * units.hertz)
    return Args(zeta=zeta, omega=omega)


def test_law(test_args: Args) -> None:
    result = decay_law.calculate_exponential_decay_constant(test_args.omega, test_args.zeta)
    assert_equal(result, 3.0 / units.second)


def test_bad_frequency(test_args: Args) -> None:
    omega_bad = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        decay_law.calculate_exponential_decay_constant(omega_bad, test_args.zeta)
    with raises(TypeError):
        decay_law.calculate_exponential_decay_constant(100, test_args.zeta)


def test_bad_damping_ratio(test_args: Args) -> None:
    zeta_bad = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        decay_law.calculate_exponential_decay_constant(test_args.omega, zeta_bad)
