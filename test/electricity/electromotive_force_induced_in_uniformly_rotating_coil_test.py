from collections import namedtuple
from pytest import fixture, raises
from sympy import pi
from symplyphysics import (assert_equal, units, Quantity, errors)
from symplyphysics.laws.electricity import electromotive_force_induced_in_uniformly_rotating_coil as electromotive_force_law

# Description
## In a homogeneous magnetostatic field with an induction equal to 0.1 tesla, a frame containing 1000 turns
## rotates uniformly with a frequency of 2 * pi * 10 radian per second. The area of the frame is 150 [centimeter^2],
## the time is 8.3e-3 second. Then the voltage value will be -46.95 volt.
## https://remote.misis.ru/courses/168/pages/13-dot-7-primiery-rieshieniia-zadach

Args = namedtuple("Args",
    ["number_turns", "induction", "contour_area", "rotation_frequency", "time"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    number_turns = 1000
    induction = Quantity(0.1 * units.tesla)
    contour_area = Quantity(150 * units.centimeter**2)
    rotation_frequency = Quantity(2 * pi * 10 * units.radian / units.second)
    time = Quantity(8.3e-3 * units.second)
    return Args(number_turns=number_turns,
        induction=induction,
        contour_area=contour_area,
        rotation_frequency=rotation_frequency,
        time=time)


def test_basic_voltage(test_args: Args) -> None:
    result = electromotive_force_law.calculate_voltage(test_args.number_turns, test_args.induction,
        test_args.contour_area, test_args.rotation_frequency, test_args.time)
    assert_equal(result, -46.95 * units.volt)


def test_bad_number_turns(test_args: Args) -> None:
    number_turns = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        electromotive_force_law.calculate_voltage(number_turns, test_args.induction,
            test_args.contour_area, test_args.rotation_frequency, test_args.time)


def test_bad_induction(test_args: Args) -> None:
    induction = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        electromotive_force_law.calculate_voltage(test_args.number_turns, induction,
            test_args.contour_area, test_args.rotation_frequency, test_args.time)
    with raises(TypeError):
        electromotive_force_law.calculate_voltage(test_args.number_turns, 100,
            test_args.contour_area, test_args.rotation_frequency, test_args.time)


def test_bad_contour_area(test_args: Args) -> None:
    contour_area = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        electromotive_force_law.calculate_voltage(test_args.number_turns, test_args.induction,
            contour_area, test_args.rotation_frequency, test_args.time)
    with raises(TypeError):
        electromotive_force_law.calculate_voltage(test_args.number_turns, test_args.induction, 100,
            test_args.rotation_frequency, test_args.time)


def test_bad_rotation_frequency(test_args: Args) -> None:
    rotation_frequency = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        electromotive_force_law.calculate_voltage(test_args.number_turns, test_args.induction,
            test_args.contour_area, rotation_frequency, test_args.time)
    with raises(TypeError):
        electromotive_force_law.calculate_voltage(test_args.number_turns, test_args.induction,
            test_args.contour_area, 100, test_args.time)


def test_bad_time(test_args: Args) -> None:
    time = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        electromotive_force_law.calculate_voltage(test_args.number_turns, test_args.induction,
            test_args.contour_area, test_args.rotation_frequency, time)
    with raises(TypeError):
        electromotive_force_law.calculate_voltage(test_args.number_turns, test_args.induction,
            test_args.contour_area, test_args.rotation_frequency, 100)
