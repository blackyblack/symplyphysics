from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.kinematics.damped_oscillations import (damped_angular_frequency as
    damped_law)

# Description
## A damped oscillator has a damping ratio of 0.3 and its undamped angular frequency
## is 100 Hz. Then its damped angular frequency amounts to 95.4 Hz.

Args = namedtuple("Args", "omega zeta")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    omega = Quantity(100.0 * units.hertz)
    zeta = 0.3
    return Args(omega=omega, zeta=zeta)


def test_law(test_args: Args) -> None:
    result = damped_law.calculate_damped_angular_frequency(test_args.omega, test_args.zeta)
    assert_equal(result, 95.4 * units.hertz)


def test_bad_frequency(test_args: Args) -> None:
    omega_bad = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        damped_law.calculate_damped_angular_frequency(omega_bad, test_args.zeta)
    with raises(TypeError):
        damped_law.calculate_damped_angular_frequency(100, test_args.zeta)


def test_bad_damping_ratio(test_args: Args) -> None:
    zeta_bad = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        damped_law.calculate_damped_angular_frequency(test_args.omega, zeta_bad)
