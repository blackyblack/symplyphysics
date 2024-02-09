from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import assert_equal, Quantity, units, errors
from symplyphysics.laws.thermodynamics import prandtl_number

# Example from https://www.omnicalculator.com/physics/prandtl-number

Args = namedtuple("Args", ["heat_capacity", "dynamic_viscosity", "thermal_conductivity"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    heat_capacity = Quantity(1.006 * units.joule / (units.kilogram * units.kelvin))
    dynamic_viscosity = Quantity(0.0000182 * units.pascal * units.second)
    thermal_conductivity = Quantity(0.025596 * units.watt / units.meter / units.kelvin)
    return Args(heat_capacity=heat_capacity,
        dynamic_viscosity=dynamic_viscosity,
        thermal_conductivity=thermal_conductivity)


def test_basic_prandtl_number(test_args: Args) -> None:
    result = prandtl_number.calculate_prandtl_number(test_args.heat_capacity,
        test_args.dynamic_viscosity, test_args.thermal_conductivity)
    assert_equal(result, 0.000715)


def test_bad_heat_capacity(test_args: Args) -> None:
    bc = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        prandtl_number.calculate_prandtl_number(
            bc,
            test_args.dynamic_viscosity,
            test_args.thermal_conductivity,
        )
    with raises(TypeError):
        prandtl_number.calculate_prandtl_number(
            100,
            test_args.dynamic_viscosity,
            test_args.thermal_conductivity,
        )


def test_bad_dynamic_viscosity(test_args: Args) -> None:
    bdv = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        prandtl_number.calculate_prandtl_number(
            test_args.heat_capacity,
            bdv,
            test_args.thermal_conductivity,
        )
    with raises(TypeError):
        prandtl_number.calculate_prandtl_number(
            test_args.heat_capacity,
            100,
            test_args.thermal_conductivity,
        )


def test_bad_thermal_conductivity(test_args: Args) -> None:
    btk = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        prandtl_number.calculate_prandtl_number(
            test_args.heat_capacity,
            test_args.dynamic_viscosity,
            btk,
        )
    with raises(TypeError):
        prandtl_number.calculate_prandtl_number(
            test_args.heat_capacity,
            test_args.dynamic_viscosity,
            100,
        )
