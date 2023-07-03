from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    errors,
    units,
    Quantity,
    SI,
    convert_to,
)
from symplyphysics.laws.thermodynamics import pressure_from_temperature_and_volume as ideal_gas_law


@fixture(name="test_args")
def test_args_fixture():
    V = Quantity(22.414 * units.liter)
    t = Quantity(273.15 * units.kelvin)
    n = Quantity(1 * units.mole)
    Args = namedtuple("Args", ["V", "t", "n"])
    return Args(V=V, t=t, n=n)


# The volume of 1 mol of any gas at STP (Standard temperature, 273.15 K and pressure, 1 atm) is measured to be 22.414 L.
def test_basic_pressure(test_args):
    result = ideal_gas_law.calculate_pressure(test_args.V, test_args.t, test_args.n)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.pressure)
    result_pressure = convert_to(result, units.pascal).evalf(2)
    assert result_pressure == approx(101.325e+3, 0.01)
    # Also check that calculated pressure = 1 atmosphere
    result_pressure_atms = convert_to(result, units.atm).evalf(2)
    assert result_pressure_atms == approx(1.0, 0.01)


def test_bad_volume(test_args):
    Vb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        ideal_gas_law.calculate_pressure(Vb, test_args.t, test_args.n)
    with raises(TypeError):
        ideal_gas_law.calculate_pressure(100, test_args.t, test_args.n)


def test_bad_temperature(test_args):
    tb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        ideal_gas_law.calculate_pressure(test_args.V, tb, test_args.n)
    with raises(TypeError):
        ideal_gas_law.calculate_pressure(test_args.V, 100, test_args.n)


def test_bad_mole_count(test_args):
    nb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        ideal_gas_law.calculate_pressure(test_args.V, test_args.t, nb)
    with raises(TypeError):
        ideal_gas_law.calculate_pressure(test_args.V, test_args.t, 100)
