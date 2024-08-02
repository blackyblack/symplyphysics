from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.thermodynamics import (
    internal_energy_change_of_ideal_gas_via_temperature as energy_law,)

# Description
## The internal energy change of a gas with isochoric heat capacity C_V = 100 J/K and a temperature
## change of 0.001 K is U = 0.1 J.

Args = namedtuple("Args", "cv dt")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    cv = Quantity(100 * units.joule / units.kelvin)
    dt = Quantity(0.001 * units.kelvin)
    return Args(cv=cv, dt=dt)


def test_law(test_args: Args) -> None:
    result = energy_law.calculate_internal_energy(test_args.cv, test_args.dt)
    assert_equal(result, 0.1 * units.joule)


def test_bad_heat_capacity(test_args: Args) -> None:
    cb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        energy_law.calculate_internal_energy(cb, test_args.dt)
    with raises(TypeError):
        energy_law.calculate_internal_energy(100, test_args.dt)


def test_bad_temperature(test_args: Args) -> None:
    tb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        energy_law.calculate_internal_energy(test_args.cv, tb)
    with raises(TypeError):
        energy_law.calculate_internal_energy(test_args.cv, 100)
