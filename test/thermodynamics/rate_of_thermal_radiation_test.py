from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.thermodynamics import rate_of_thermal_radiation as radiation_law

# Description
## The rate of thermal emission of a body with surface emissivity epsilon = 0.5, area A = 1 m**2
## and temperature T = 300 K is approximately 229 W.

Args = namedtuple("Args", "e a t")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    e = 0.5
    a = Quantity(1 * units.meter**2)
    t = Quantity(300 * units.kelvin)
    return Args(e=e, a=a, t=t)


def test_law(test_args: Args) -> None:
    result = radiation_law.calculate_energy_radiation_rate(test_args.e, test_args.a, test_args.t)
    assert_equal(result, 229 * units.watt, tolerance=3e-3)


def test_bad_emissivity(test_args: Args) -> None:
    eb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        radiation_law.calculate_energy_radiation_rate(eb, test_args.a, test_args.t)


def test_bad_area(test_args: Args) -> None:
    ab = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        radiation_law.calculate_energy_radiation_rate(test_args.e, ab, test_args.t)
    with raises(TypeError):
        radiation_law.calculate_energy_radiation_rate(test_args.e, 100, test_args.t)


def test_bad_temperature(test_args: Args) -> None:
    tb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        radiation_law.calculate_energy_radiation_rate(test_args.e, test_args.a, tb)
    with raises(TypeError):
        radiation_law.calculate_energy_radiation_rate(test_args.e, test_args.a, 100)
