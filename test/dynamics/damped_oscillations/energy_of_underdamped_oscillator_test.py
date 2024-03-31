from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.dynamics.damped_oscillations import (energy_of_underdamped_oscillator as
    energy_law)

# Description
## An object is performing damped oscillating motion. It has mass m = 50 g, the undamped angular
## frequency of the oscillator is w0 = 5 rad/s, the maximum amplitude of the damped oscillations
## is A = 0.3 m, and the exponential decay constant of the oscillator is lambda = 0.3 rad/s.
## At t = 1 s, the energy of the oscillator E ~= 0.0309 J.

Args = namedtuple("Args", "m a w0 lambda_ t")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    m = Quantity(50.0 * units.gram)
    a = Quantity(0.3 * units.meter)
    w0 = Quantity(5.0 * units.radian / units.second)
    lambda_ = Quantity(0.3 * units.radian / units.second)
    t = Quantity(1.0 * units.second)
    return Args(m=m, a=a, w0=w0, lambda_=lambda_, t=t)


def test_law(test_args: Args) -> None:
    result = energy_law.calculate_oscillator_energy(test_args.m, test_args.a, test_args.w0,
        test_args.lambda_, test_args.t)
    assert_equal(result, 0.0309 * units.joule)


def test_bad_mass(test_args: Args) -> None:
    mb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        energy_law.calculate_oscillator_energy(mb, test_args.a, test_args.w0, test_args.lambda_,
            test_args.t)
    with raises(TypeError):
        energy_law.calculate_oscillator_energy(100, test_args.a, test_args.w0, test_args.lambda_,
            test_args.t)


def test_bad_amplitude(test_args: Args) -> None:
    ab = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        energy_law.calculate_oscillator_energy(test_args.m, ab, test_args.w0, test_args.lambda_,
            test_args.t)
    with raises(TypeError):
        energy_law.calculate_oscillator_energy(test_args.m, 100, test_args.w0, test_args.lambda_,
            test_args.t)


def test_bad_frequency(test_args: Args) -> None:
    wb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        energy_law.calculate_oscillator_energy(test_args.m, test_args.a, wb, test_args.lambda_,
            test_args.t)
    with raises(TypeError):
        energy_law.calculate_oscillator_energy(test_args.m, test_args.a, 100, test_args.lambda_,
            test_args.t)


def test_bad_decay_constant(test_args: Args) -> None:
    lb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        energy_law.calculate_oscillator_energy(test_args.m, test_args.a, test_args.w0, lb,
            test_args.t)
    with raises(TypeError):
        energy_law.calculate_oscillator_energy(test_args.m, test_args.a, test_args.w0, 100,
            test_args.t)


def test_bad_time(test_args: Args) -> None:
    tb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        energy_law.calculate_oscillator_energy(test_args.m, test_args.a, test_args.w0,
            test_args.lambda_, tb)
    with raises(TypeError):
        energy_law.calculate_oscillator_energy(test_args.m, test_args.a, test_args.w0,
            test_args.lambda_, 100)
