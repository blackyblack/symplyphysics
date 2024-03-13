from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    prefixes,
    Quantity,
)
from symplyphysics.laws.waves import pressure_amplitude_in_sound_wave as pressure_law

# Description
## A sound wave of angular frequency w = 6.28e3 rad/s is propagating in air of density
## rho = 1.21 kg/m**3 at a speed of 343 m/s. The maximum displacement of the air particles
## is 11 µm. Thus the amplitude of the pressure change is 28 Pa.

Args = namedtuple("Args", "v rho w a")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    v = Quantity(343 * units.meter / units.second)
    rho = Quantity(1.21 * units.kilogram / units.meter**3)
    w = Quantity(6.23e3 * units.radian / units.second)
    a = Quantity(11 * prefixes.micro * units.meter)
    return Args(v=v, rho=rho, w=w, a=a)


def test_law(test_args: Args) -> None:
    result = pressure_law.calculate_pressure_amplitude(test_args.v, test_args.rho, test_args.w, test_args.a)
    assert_equal(result, 28 * units.pascal, tolerance=2e-2)


def test_bad_speed(test_args: Args) -> None:
    vb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        pressure_law.calculate_pressure_amplitude(vb, test_args.rho, test_args.w, test_args.a)
    with raises(TypeError):
        pressure_law.calculate_pressure_amplitude(100, test_args.rho, test_args.w, test_args.a)


def test_bad_density(test_args: Args) -> None:
    rho_bad = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        pressure_law.calculate_pressure_amplitude(test_args.v, rho_bad, test_args.w, test_args.a)
    with raises(TypeError):
        pressure_law.calculate_pressure_amplitude(test_args.v, 100, test_args.w, test_args.a)


def test_bad_angular_frequency(test_args: Args) -> None:
    wb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        pressure_law.calculate_pressure_amplitude(test_args.v, test_args.rho, wb, test_args.a)
    with raises(TypeError):
        pressure_law.calculate_pressure_amplitude(test_args.v, test_args.rho, 100, test_args.a)


def test_bad_displacement_amplitude(test_args: Args) -> None:
    ab = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        pressure_law.calculate_pressure_amplitude(test_args.v, test_args.rho, test_args.w, ab)
    with raises(TypeError):
        pressure_law.calculate_pressure_amplitude(test_args.v, test_args.rho, test_args.w, 100)
