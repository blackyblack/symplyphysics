from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    errors,
    units,
    Quantity,
    assert_equal
)

from symplyphysics.laws.hydro import flow_of_liquid_from_viscosity_and_pressure as flow_of_liquid

## Source of numbers: https://testbook.com/physics-formulas/poiseuilles-law-formula

Args = namedtuple("Args", ["pressure_differential", "radius", "viscosity", "length",])


@fixture(name="test_args")
def test_args_fixture():
    pressure_differential = Quantity(0.14 * units.pascal)
    radius = Quantity(3 * units.meter)
    viscosity = Quantity(0.056 * units.pascal * units.second)
    length = Quantity(8 * units.meter)
    return Args(
        pressure_differential=pressure_differential,
        radius=radius,
        viscosity=viscosity,
        length=length,
    )


def test_basic_law(test_args: Args) -> None:
    result = flow_of_liquid.calculate_flow_of_liquid(test_args.pressure_differential, test_args.radius, test_args.viscosity, test_args.length)
    assert_equal(result, 9.94 * (units.meter**3 / units.second))


def test_bad_pressure(test_args: Args) -> None:
    bad_pressure = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        flow_of_liquid.calculate_flow_of_liquid(bad_pressure, test_args.radius, test_args.viscosity, test_args.length)
    with raises(TypeError):
        flow_of_liquid.calculate_flow_of_liquid(100, test_args.radius, test_args.viscosity, test_args.length)


def test_bad_length(test_args: Args) -> None:
    bad_length = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        flow_of_liquid.calculate_flow_of_liquid(test_args.pressure_differential, bad_length, test_args.viscosity, test_args.length)
    with raises(TypeError):
        flow_of_liquid.calculate_flow_of_liquid(test_args.pressure_differential, 100, test_args.viscosity, test_args.length)
    with raises(errors.UnitsError):
        flow_of_liquid.calculate_flow_of_liquid(test_args.pressure_differential, test_args.radius, test_args.viscosity, bad_length)
    with raises(TypeError):
        flow_of_liquid.calculate_flow_of_liquid(test_args.pressure_differential, test_args.radius, test_args.viscosity, 100)


def test_bad_viscosity(test_args: Args) -> None:
    bad_viscosity = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        flow_of_liquid.calculate_flow_of_liquid(test_args.pressure_differential, test_args.radius, bad_viscosity, test_args.length)
    with raises(TypeError):
        flow_of_liquid.calculate_flow_of_liquid(test_args.pressure_differential, test_args.radius, 100, test_args.length)
