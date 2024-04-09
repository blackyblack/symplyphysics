from collections import namedtuple
from pytest import fixture, raises
from sympy import sqrt
from symplyphysics import (errors, units, Quantity, assert_equal)

from symplyphysics.laws.electricity.circuits import scattering_matrix as parameters_matrix_law

## The reflected wave coefficient per input port is 0.2 [watt^(1/2)], the reflected wave coefficient per output port is 0.2 [watt^(1/2)].
## S-parameters S11, S12, S21, S22 are equal to 0.2, 0.5, 0.5, 0.2, respectively.
## Then the values of the incident wave coefficient per input port and incident wave coefficient per output port will be equal
## to 0.286 [watt^(1/2)].

Args = namedtuple("Args", ["input_reflected_wave", "output_reflected_wave", "parameters"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    input_reflected_wave = Quantity(0.2 * sqrt(units.watt))
    output_reflected_wave = Quantity(0.2 * sqrt(units.watt))
    parameters = ((0.2, 0.5), (0.5, 0.2))
    return Args(input_reflected_wave=input_reflected_wave,
        output_reflected_wave=output_reflected_wave,
        parameters=parameters
        )


def test_basic_waves(test_args: Args) -> None:
    result = parameters_matrix_law.calculate_waves(test_args.input_reflected_wave, test_args.output_reflected_wave, test_args.parameters)
    assert_equal(result[0], 0.286 * sqrt(units.watt))
    assert_equal(result[1], 0.286 * sqrt(units.watt))


def test_bad_waves(test_args: Args) -> None:
    bad_wave = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        parameters_matrix_law.calculate_waves(bad_wave, test_args.output_reflected_wave, test_args.parameters)
    with raises(TypeError):
        parameters_matrix_law.calculate_waves(100, test_args.output_reflected_wave, test_args.parameters)
    with raises(errors.UnitsError):
        parameters_matrix_law.calculate_waves(test_args.input_reflected_wave, bad_wave, test_args.parameters)
    with raises(TypeError):
        parameters_matrix_law.calculate_waves(test_args.input_reflected_wave, 100, test_args.parameters)
