from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.waves import (
    phase_velocity_from_angular_frequency_and_wavenumber as phase_velocity_law,)

# Description
## A wave has angular frequency w = 9 Hz and angular wavenumber k = 0.3 1/m.
## Its phase speed is 30 m/s.

Args = namedtuple("Args", "w k")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    w = Quantity(9.0 * units.hertz)
    k = Quantity(0.3 / units.meter)
    return Args(w=w, k=k)


def test_law(test_args: Args) -> None:
    result = phase_velocity_law.calculate_phase_velocity(test_args.w, test_args.k)
    assert_equal(result, 30.0 * units.meter / units.second)


def test_bad_frequency(test_args: Args) -> None:
    wb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        phase_velocity_law.calculate_phase_velocity(wb, test_args.k)
    with raises(TypeError):
        phase_velocity_law.calculate_phase_velocity(100, test_args.k)


def test_bad_wavenumber(test_args: Args) -> None:
    kb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        phase_velocity_law.calculate_phase_velocity(test_args.w, kb)
    with raises(TypeError):
        phase_velocity_law.calculate_phase_velocity(test_args.w, 100)
