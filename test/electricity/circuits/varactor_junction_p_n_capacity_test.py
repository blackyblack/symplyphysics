from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, units, Quantity, errors, prefixes)
from symplyphysics.laws.electricity.circuits import varactor_junction_p_n_capacity as capacity_law

# Description
## The capacity of the varactor without bias voltage is equal to 10 nanofarad. The voltage is 0.5 volt,
## the material parameter is 1 volt. The doping coefficient is 0.5. Then the capacity of the varactor
## at a given voltage will be equal to 14.14 nanofarad.

Args = namedtuple("Args", ["junction_capacitance_without_bias_voltage", "voltage", "material_parameter", "doping_coefficient"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    junction_capacitance_without_bias_voltage = Quantity(10 * prefixes.nano * units.farad)
    voltage = Quantity(0.5 * units.volt)
    material_parameter = Quantity(1 * units.volt)
    doping_coefficient = 0.5

    return Args(junction_capacitance_without_bias_voltage=junction_capacitance_without_bias_voltage,
                voltage=voltage,
                material_parameter=material_parameter,
                doping_coefficient=doping_coefficient)


def test_basic_junction_capacitance(test_args: Args) -> None:
    result = capacity_law.calculate_junction_capacitance(test_args.junction_capacitance_without_bias_voltage, test_args.voltage, test_args.material_parameter, test_args.doping_coefficient)
    assert_equal(result, 14.14 * prefixes.nano * units.farad)


def test_bad_junction_capacitance_without_bias_voltage(test_args: Args) -> None:
    junction_capacitance_without_bias_voltage = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        capacity_law.calculate_junction_capacitance(junction_capacitance_without_bias_voltage, test_args.voltage, test_args.material_parameter, test_args.doping_coefficient)
    with raises(TypeError):
        capacity_law.calculate_junction_capacitance(100, test_args.voltage, test_args.material_parameter, test_args.doping_coefficient)


def test_bad_voltage(test_args: Args) -> None:
    voltage = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        capacity_law.calculate_junction_capacitance(test_args.junction_capacitance_without_bias_voltage, voltage, test_args.material_parameter, test_args.doping_coefficient)
    with raises(TypeError):
        capacity_law.calculate_junction_capacitance(test_args.junction_capacitance_without_bias_voltage, 100, test_args.material_parameter, test_args.doping_coefficient)


def test_bad_material_parameter(test_args: Args) -> None:
    material_parameter = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        capacity_law.calculate_junction_capacitance(test_args.junction_capacitance_without_bias_voltage, test_args.voltage, material_parameter, test_args.doping_coefficient)
    with raises(TypeError):
        capacity_law.calculate_junction_capacitance(test_args.junction_capacitance_without_bias_voltage, test_args.voltage, 100, test_args.doping_coefficient)


def test_bad_doping_coefficient(test_args: Args) -> None:
    doping_coefficient = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        capacity_law.calculate_junction_capacitance(test_args.junction_capacitance_without_bias_voltage, test_args.voltage, test_args.material_parameter, doping_coefficient)
