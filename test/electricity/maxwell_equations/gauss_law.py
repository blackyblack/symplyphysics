from collections import namedtuple
from pytest import fixture, raises
from sympy import sin
from symplyphysics import (
    errors,
    units,
    Quantity,
    assert_equal,
)
from symplyphysics.core.coordinate_systems.coordinate_systems import CoordinateSystem
from symplyphysics.core.fields.vector_field import VectorField
from sympy.physics.units import electric_constant

from symplyphysics.laws.electricity.maxwell_equations import gauss_law

## TODO

Args = namedtuple("Args", ["electrical_induction", "x", "y", "z"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    C = CoordinateSystem()
    electric_intensity = Quantity(1e4 * units.volt / units.meter)
    x = Quantity(0.493 * units.meter)
    y = Quantity(0.5 * units.meter)
    z = Quantity(0.01249 * units.meter)
    electrical_induction = VectorField([sin(x / Quantity(1 * units.meter)) * electric_intensity * electric_constant,
                    sin(y / Quantity(1 * units.meter)) * electric_intensity * electric_constant,
                    sin(z / Quantity(1 * units.meter)) * electric_intensity * electric_constant], C)
    return Args(electrical_induction=electrical_induction,
        x=x,
        y=y,
        z=z)


def test_basic_bulk_density_of_electric_charge(test_args: Args) -> None:
    result = gauss_law.calculate_bulk_density_of_electric_charge(test_args.electrical_induction, test_args.x, test_args.y, test_args.z)
    assert_equal(result, 1 * units.coulomb / units.meter**3)


def test_bad_coordinates(test_args: Args) -> None:
    bad_coordinate = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        gauss_law.calculate_bulk_density_of_electric_charge(test_args.electrical_induction, bad_coordinate, test_args.y, test_args.z)
    with raises(errors.UnitsError):
        gauss_law.calculate_bulk_density_of_electric_charge(test_args.electrical_induction, test_args.x, bad_coordinate, test_args.z)
    with raises(errors.UnitsError):
        gauss_law.calculate_bulk_density_of_electric_charge(test_args.electrical_induction, test_args.x, test_args.y, bad_coordinate)
    with raises(TypeError):
        gauss_law.calculate_bulk_density_of_electric_charge(test_args.electrical_induction, 100, test_args.y, test_args.z)
    with raises(TypeError):
        gauss_law.calculate_bulk_density_of_electric_charge(test_args.electrical_induction, test_args.x, 100, test_args.z)
    with raises(TypeError):
        gauss_law.calculate_bulk_density_of_electric_charge(test_args.electrical_induction, test_args.x, test_args.y, 100)
