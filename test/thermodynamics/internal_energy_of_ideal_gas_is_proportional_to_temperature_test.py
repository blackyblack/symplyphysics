from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    prefixes,
    units,
    Quantity,
)
from symplyphysics.laws.thermodynamics import (
    internal_energy_of_ideal_gas_is_proportional_to_temperature as energy_law,
)

# Description
## The internal energy of a gas with isochoric heat capacity C_V = 100 J/K and a temperature
## of 100 K is U = 10 kJ.

Args = namedtuple("Args", "cv t")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    cv = Quantity(100 * units.joule / units.kelvin)
    t = Quantity(100 * units.kelvin)
    return Args(cv=cv, t=t)


def test_law(test_args: Args) -> None:
    result = energy_law.calculate_internal_energy(test_args.cv, test_args.t)
    assert_equal(result, 10 * prefixes.kilo * units.joule)


def test_bad_heat_capacity(test_args: Args) -> None:
    cb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        energy_law.calculate_internal_energy(cb, test_args.t)
    with raises(TypeError):
        energy_law.calculate_internal_energy(100, test_args.t)


def test_bad_temperature(test_args: Args) -> None:
    tb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        energy_law.calculate_internal_energy(test_args.cv, tb)
    with raises(TypeError):
        energy_law.calculate_internal_energy(test_args.cv, 100)