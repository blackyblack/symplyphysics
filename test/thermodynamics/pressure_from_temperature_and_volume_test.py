from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.thermodynamics import pressure_from_temperature_and_volume as ideal_gas_law

Args = namedtuple("Args", ["V", "t", "n"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    V = Quantity(22.414 * units.liter)
    t = Quantity(273.15 * units.kelvin)
    n = Quantity(1 * units.mole)
    return Args(V=V, t=t, n=n)


# The volume of 1 mol of any gas at STP (Standard temperature, 273.15 K and pressure, 1 atm) is measured to be 22.414 L.
def test_basic_pressure(test_args: Args) -> None:
    result = ideal_gas_law.calculate_pressure(test_args.V, test_args.t, test_args.n)
    assert_equal(result, 1.01325e5 * units.pascal)
    # Also check that calculated pressure = 1 atmosphere
    assert_equal(result, 1 * units.atm)


def test_bad_volume(test_args: Args) -> None:
    Vb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        ideal_gas_law.calculate_pressure(Vb, test_args.t, test_args.n)
    with raises(TypeError):
        ideal_gas_law.calculate_pressure(100, test_args.t, test_args.n)


def test_bad_temperature(test_args: Args) -> None:
    tb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        ideal_gas_law.calculate_pressure(test_args.V, tb, test_args.n)
    with raises(TypeError):
        ideal_gas_law.calculate_pressure(test_args.V, 100, test_args.n)


def test_bad_mole_count(test_args: Args) -> None:
    nb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        ideal_gas_law.calculate_pressure(test_args.V, test_args.t, nb)
    with raises(TypeError):
        ideal_gas_law.calculate_pressure(test_args.V, test_args.t, 100)
