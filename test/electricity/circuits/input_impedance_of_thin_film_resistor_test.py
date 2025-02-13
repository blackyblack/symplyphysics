from collections import namedtuple
from pytest import fixture, raises
from sympy import I, pi
from symplyphysics import (assert_equal, units, Quantity, errors, prefixes)
from symplyphysics.laws.electricity.circuits import input_impedance_of_thin_film_resistor as impedance_law

# Description
## The resistance and capacitance of a thin-film resistor are 40 ohm and 40 nanofarad, respectively. The angular frequency is 2*pi*10^6 radian/second.
## Then the input impedance of the thin-film resistor is (3.27 - 10.96 * I) ohm.

Args = namedtuple("Args", [
    "resistance",
    "angular_frequency",
    "capacitance",
])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    resistance = Quantity(40 * units.ohm)
    angular_frequency = Quantity(2 * pi * prefixes.mega * units.hertz)
    capacitance = Quantity(40 * prefixes.nano * units.farad)

    return Args(
        resistance=resistance,
        angular_frequency=angular_frequency,
        capacitance=capacitance,
    )


def test_basic_input_impedance(test_args: Args) -> None:
    result = impedance_law.calculate_input_impedance(test_args.resistance, test_args.angular_frequency,
        test_args.capacitance)
    assert_equal(result, (3.27 - 10.96 * I) * units.ohm)


def test_bad_resistance(test_args: Args) -> None:
    resistance = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        impedance_law.calculate_input_impedance(resistance, test_args.angular_frequency,
            test_args.capacitance)
    with raises(TypeError):
        impedance_law.calculate_input_impedance(100, test_args.angular_frequency, test_args.capacitance)


def test_bad_angular_frequency(test_args: Args) -> None:
    angular_frequency = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        impedance_law.calculate_input_impedance(test_args.resistance, angular_frequency,
            test_args.capacitance)
    with raises(TypeError):
        impedance_law.calculate_input_impedance(test_args.resistance, 100, test_args.capacitance)


def test_bad_capacitance(test_args: Args) -> None:
    capacitance = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        impedance_law.calculate_input_impedance(test_args.resistance, test_args.angular_frequency,
            capacitance)
    with raises(TypeError):
        impedance_law.calculate_input_impedance(test_args.resistance, test_args.angular_frequency, 100)
