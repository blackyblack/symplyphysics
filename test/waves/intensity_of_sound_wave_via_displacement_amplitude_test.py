from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    prefixes,
    Quantity,
)
from symplyphysics.laws.waves import intensity_of_sound_wave_via_displacement_amplitude as intensity_law

# Description
## A sound wave is propagating in air (rho = 1.21 kg/m**3) with a speed of 343 m/s at an angular
## frequency of 6.28e3 rad/s. The displacement amplitude of the sound wave is 10 µm. The intensity
## of the wave is 0.82 W/m**2.

Args = namedtuple("Args", "rho v w s_max")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    rho = Quantity(1.21 * units.kilogram / units.meter**3)
    v = Quantity(343 * units.meter / units.second)
    w = Quantity(6.28e3 * units.radian / units.second)
    s_max = Quantity(10 * prefixes.micro * units.meter)
    return Args(rho=rho, v=v, w=w, s_max=s_max)


def test_law(test_args: Args) -> None:
    result = intensity_law.calculate_intensity(test_args.rho, test_args.v, test_args.w,
        test_args.s_max)
    assert_equal(result, 0.82 * units.watt / units.meter**2, tolerance=2e-3)


def test_bad_density(test_args: Args) -> None:
    rho_bad = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        intensity_law.calculate_intensity(rho_bad, test_args.v, test_args.w, test_args.s_max)
    with raises(TypeError):
        intensity_law.calculate_intensity(100, test_args.v, test_args.w, test_args.s_max)


def test_bad_speed(test_args: Args) -> None:
    vb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        intensity_law.calculate_intensity(test_args.rho, vb, test_args.w, test_args.s_max)
    with raises(TypeError):
        intensity_law.calculate_intensity(test_args.rho, 100, test_args.w, test_args.s_max)


def test_bad_angular_frequency(test_args: Args) -> None:
    wb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        intensity_law.calculate_intensity(test_args.rho, test_args.v, wb, test_args.s_max)
    with raises(TypeError):
        intensity_law.calculate_intensity(test_args.rho, test_args.v, 100, test_args.s_max)


def test_bad_displacement(test_args: Args) -> None:
    sb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        intensity_law.calculate_intensity(test_args.rho, test_args.v, test_args.w, sb)
    with raises(TypeError):
        intensity_law.calculate_intensity(test_args.rho, test_args.v, test_args.w, 100)
