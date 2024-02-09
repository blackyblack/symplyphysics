from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.electricity import amount_energy_from_voltage_time_resistance as joule_lenz_law

#  How much energy can a household electric kettle produce to heat water in one minute
#  at 220 volts and a heater resistance of 36 ohms?

Args = namedtuple("Args", ["U", "t", "R"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    U = Quantity(220 * units.volt)
    t = Quantity(60 * units.second)
    R = Quantity(36 * units.ohm)
    return Args(U=U, t=t, R=R)


def test_basic_amount(test_args: Args) -> None:
    result = joule_lenz_law.calculate_amount_energy(test_args.U, test_args.t, test_args.R)
    assert_equal(result, 80666.6 * units.joule)


def test_bad_voltage(test_args: Args) -> None:
    Ub = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        joule_lenz_law.calculate_amount_energy(Ub, test_args.t, test_args.R)
    with raises(TypeError):
        joule_lenz_law.calculate_amount_energy(100, test_args.t, test_args.R)


def test_bad_time(test_args: Args) -> None:
    tb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        joule_lenz_law.calculate_amount_energy(test_args.U, tb, test_args.R)
    with raises(TypeError):
        joule_lenz_law.calculate_amount_energy(test_args.U, 100, test_args.R)


def test_bad_resistance(test_args: Args) -> None:
    Rb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        joule_lenz_law.calculate_amount_energy(test_args.U, test_args.t, Rb)
    with raises(TypeError):
        joule_lenz_law.calculate_amount_energy(test_args.U, test_args.t, 100)
