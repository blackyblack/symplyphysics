from collections import namedtuple
from pytest import fixture, raises
from sympy import pi
from symplyphysics import (assert_equal, units, Quantity, errors)
from symplyphysics.laws.electricity import electromotive_force_induced_in_rotating_rod as voltage_law

# Description
## The length of the rod is 0.1 meter, the magnetic field induction is 0.4 tesla, and the rotation
## frequency of the rod is 2 * pi * 16 radian per second. Then the voltage of the rod will be 0.201 volt.
## https://www.chertov.org.ua/viev_zadachi.php?num=11&par=25

Args = namedtuple("Args", ["magnetic_induction", "rotation_frequency", "rod_length"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    magnetic_induction = Quantity(0.4 * units.tesla)
    rotation_frequency = Quantity(2 * pi * 16 * (units.radian / units.second))
    rod_length = Quantity(0.1 * units.meter)

    return Args(magnetic_induction=magnetic_induction,
        rotation_frequency=rotation_frequency,
        rod_length=rod_length)


def test_basic_voltage(test_args: Args) -> None:
    result = voltage_law.calculate_voltage(test_args.magnetic_induction,
        test_args.rotation_frequency, test_args.rod_length)
    assert_equal(result, 0.201 * units.volt)


def test_bad_magnetic_induction(test_args: Args) -> None:
    magnetic_induction = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        voltage_law.calculate_voltage(magnetic_induction, test_args.rotation_frequency,
            test_args.rod_length)
    with raises(TypeError):
        voltage_law.calculate_voltage(100, test_args.rotation_frequency, test_args.rod_length)


def test_bad_rotation_frequency(test_args: Args) -> None:
    rotation_frequency = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        voltage_law.calculate_voltage(test_args.magnetic_induction, rotation_frequency,
            test_args.rod_length)
    with raises(TypeError):
        voltage_law.calculate_voltage(test_args.magnetic_induction, 100, test_args.rod_length)


def test_bad_rod_length(test_args: Args) -> None:
    rod_length = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        voltage_law.calculate_voltage(test_args.magnetic_induction, test_args.rotation_frequency,
            rod_length)
    with raises(TypeError):
        voltage_law.calculate_voltage(test_args.magnetic_induction, test_args.rotation_frequency,
            100)
