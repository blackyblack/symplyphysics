from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)

from symplyphysics.laws.thermodynamics import dynamic_viscosity_of_gas_from_control_viscosity_and_temperature as dynamic_viscosity
## Link of numbers: https://ru.wikipedia.org/wiki/Вязкость#Вязкость_воды


@fixture(name="test_args")
def test_args_fixture():
    set_temperature = Quantity(298 * units.kelvin)
    control_temperature = Quantity(291.15 * units.kelvin)
    control_viscosity = Quantity(1.827e-5 * units.pascal * units.second)
    sutherland_constant = Quantity(120 * units.kelvin)
    Args = namedtuple("Args",
        ["set_temperature", "control_temperature", "control_viscosity", "sutherland_constant"])
    return Args(set_temperature=set_temperature,
        control_temperature=control_temperature,
        control_viscosity=control_viscosity,
        sutherland_constant=sutherland_constant)


def test_basic_law(test_args):
    result = dynamic_viscosity.calculate_dynamic_viscosity(test_args.control_viscosity,
        test_args.control_temperature, test_args.sutherland_constant, test_args.set_temperature)
    assert_equal(result, 1.86e-5 * units.pascal * units.second)


def test_bad_temperature(test_args):
    bad_temperature = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        dynamic_viscosity.calculate_dynamic_viscosity(test_args.control_viscosity, bad_temperature,
            test_args.sutherland_constant, test_args.set_temperature)
    with raises(errors.UnitsError):
        dynamic_viscosity.calculate_dynamic_viscosity(test_args.control_viscosity,
            test_args.control_temperature, bad_temperature, test_args.set_temperature)
    with raises(errors.UnitsError):
        dynamic_viscosity.calculate_dynamic_viscosity(test_args.control_viscosity,
            test_args.control_temperature, test_args.sutherland_constant, bad_temperature)
    with raises(TypeError):
        dynamic_viscosity.calculate_dynamic_viscosity(test_args.control_viscosity, 100,
            test_args.sutherland_constant, test_args.set_temperature)
    with raises(TypeError):
        dynamic_viscosity.calculate_dynamic_viscosity(test_args.control_viscosity,
            test_args.control_temperature, 100, test_args.set_temperature)
    with raises(TypeError):
        dynamic_viscosity.calculate_dynamic_viscosity(test_args.control_viscosity,
            test_args.control_temperature, test_args.sutherland_constant, 100)


def test_bad_viscosity(test_args):
    bad_viscosity = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        dynamic_viscosity.calculate_dynamic_viscosity(bad_viscosity, test_args.control_temperature,
            test_args.sutherland_constant, test_args.set_temperature)
    with raises(TypeError):
        dynamic_viscosity.calculate_dynamic_viscosity(100, test_args.control_temperature,
            test_args.sutherland_constant, test_args.set_temperature)
