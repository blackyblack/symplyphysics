from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (errors, units, Quantity, assert_equal, prefixes)

from symplyphysics.laws.electricity.transmission_lines import input_impedance_from_parameters_transmission_matrix as impedance_law

## ABCD-parameters: A is equal to 50, B is equal to 100 ohm,
## C is equal to 1000 siemens, D is equal to 50.
## Load resistance is equal to 100 ohm.
## Then the input impedance is equal to 50.97 milliohm.

Args = namedtuple("Args", ["load_resistance", "parameters"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    load_resistance = Quantity(100 * units.ohm)
    parameters = ((50, Quantity(100 * units.ohm)), (Quantity(1000 * units.siemens), 50))
    return Args(load_resistance=load_resistance, parameters=parameters)


def test_basic_input_impedance(test_args: Args) -> None:
    result = impedance_law.calculate_input_impedance(test_args.load_resistance,
        test_args.parameters)
    assert_equal(result, 50.97 * prefixes.milli * units.ohm)


def test_bad_load_resistance(test_args: Args) -> None:
    bad_load_resistance = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        impedance_law.calculate_input_impedance(bad_load_resistance, test_args.parameters)
    with raises(TypeError):
        impedance_law.calculate_input_impedance(100, test_args.parameters)


def test_bad_parameters(test_args: Args) -> None:
    bad_parameter = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        impedance_law.calculate_input_impedance(test_args.load_resistance,
            ((test_args.parameters[0][0], bad_parameter),
            (test_args.parameters[1][0], test_args.parameters[1][1])))
    with raises(TypeError):
        impedance_law.calculate_input_impedance(test_args.load_resistance,
            ((test_args.parameters[0][0], 100),
            (test_args.parameters[1][0], test_args.parameters[1][1])))
    with raises(errors.UnitsError):
        impedance_law.calculate_input_impedance(test_args.load_resistance,
            ((test_args.parameters[0][0], test_args.parameters[0][1]),
            (bad_parameter, test_args.parameters[1][1])))
    with raises(TypeError):
        impedance_law.calculate_input_impedance(test_args.load_resistance,
            ((test_args.parameters[0][0], test_args.parameters[0][1]),
            (100, test_args.parameters[1][1])))
