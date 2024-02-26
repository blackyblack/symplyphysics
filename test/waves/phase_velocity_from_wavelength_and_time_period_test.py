from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.waves import (
    phase_velocity_from_wavelength_and_time_period as phase_velocity_law,
)

# Description
## A wave is described as having a wavelength lambda = 1 m and a time period
## T = 0.5 s. Then its phase velocity is 2 m/s.

Args = namedtuple("Args", "lambda_ t")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    lambda_ = Quantity(1.0 * units.meter)
    t = Quantity(0.5 * units.second)
    return Args(lambda_=lambda_, t=t)


def test_law(test_args: Args) -> None:
    result = phase_velocity_law.calculate_phase_velocity(test_args.lambda_, test_args.t)
    assert_equal(result, 2.0 * units.meter / units.second)


def test_bad_wavelength(test_args: Args) -> None:
    lambda_bad = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        phase_velocity_law.calculate_phase_velocity(lambda_bad, test_args.t)
    with raises(TypeError):
        phase_velocity_law.calculate_phase_velocity(100, test_args.t)


def test_bad_time_period(test_args: Args) -> None:
    tb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        phase_velocity_law.calculate_phase_velocity(test_args.lambda_, tb)
    with raises(TypeError):
        phase_velocity_law.calculate_phase_velocity(test_args.lambda_, 100)
