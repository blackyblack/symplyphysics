from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.core.symbols.celsius import Celsius, to_kelvin_quantity
from symplyphysics.laws.thermodynamics.maxwell_boltzmann_statistics import (
    velocity_component_distribution as distribution_law,
)

# Description
## The value of the velocity component distribution function for an ensemble of
## Argon particles (particle mass m = 39.948 u) at equilibrium temperature T = 10°C
## at speed v_k = 10 m/s is 1.64e-3 1/(m/s).

Args = namedtuple("Args", "v m t")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    v = Quantity(10 * units.meter / units.second)
    m = Quantity(39.948 * units.amu)
    t = to_kelvin_quantity(Celsius(10))
    return Args(v=v, m=m, t=t)


def test_law(test_args: Args) -> None:
    result = distribution_law.calculate_velocity_component_distribution(test_args.v, test_args.m, test_args.t)
    assert_equal(result, 1.64e-3 / (units.meter / units.second), tolerance=2e-3)


def test_bad_velocity(test_args: Args) -> None:
    vb = Quantity(1.0 * units.second)
    with raises(errors.UnitsError):
        distribution_law.calculate_velocity_component_distribution(vb, test_args.m, test_args.t)
    with raises(TypeError):
        distribution_law.calculate_velocity_component_distribution(100, test_args.m, test_args.t)


def test_bad_mass(test_args: Args) -> None:
    mb = Quantity(1.0 * units.second)
    with raises(errors.UnitsError):
        distribution_law.calculate_velocity_component_distribution(test_args.v, mb, test_args.t)
    with raises(TypeError):
        distribution_law.calculate_velocity_component_distribution(test_args.v, 100, test_args.t)


def test_bad_temperature(test_args: Args) -> None:
    tb = Quantity(1.0 * units.second)
    with raises(errors.UnitsError):
        distribution_law.calculate_velocity_component_distribution(test_args.v, test_args.m, tb)
    with raises(TypeError):
        distribution_law.calculate_velocity_component_distribution(test_args.v, test_args.m, 100)
