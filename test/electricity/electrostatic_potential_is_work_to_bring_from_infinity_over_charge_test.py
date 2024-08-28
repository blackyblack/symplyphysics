from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, units, Quantity, errors)
from symplyphysics.laws.electricity import electrostatic_potential_is_work_to_bring_from_infinity_over_charge as work_law

# Description
## At a voltage of 5 volts and a charge of 1 coulomb, work of charge transfer
## will be equal to 5 joules.
## https://www.fxyz.ru/формулы_по_физике/электричество/цепи_постоянного_тока/работа_и_мощность_электрического_тока/работа_электрического_тока/

Args = namedtuple("Args", ["charge", "voltage"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    charge = Quantity(1 * units.coulomb)
    voltage = Quantity(5 * units.volt)
    return Args(charge=charge, voltage=voltage)


def test_basic_work_of_charge_transfer(test_args: Args) -> None:
    result = work_law.calculate_work(test_args.charge, test_args.voltage)
    assert_equal(result, 5 * units.joule)


def test_bad_charge(test_args: Args) -> None:
    charge = Quantity(1 * units.joule)
    with raises(errors.UnitsError):
        work_law.calculate_work(charge, test_args.voltage)
    with raises(TypeError):
        work_law.calculate_work(100, test_args.voltage)


def test_bad_voltage(test_args: Args) -> None:
    voltage = Quantity(1 * units.joule)
    with raises(errors.UnitsError):
        work_law.calculate_work(test_args.charge, voltage)
    with raises(TypeError):
        work_law.calculate_work(test_args.charge, 100)
