from collections import namedtuple
from pytest import fixture, raises
from sympy import sin
from symplyphysics import units, Quantity, assert_equal, errors, quantities
from symplyphysics.laws.electricity.maxwell_equations import (
    charge_density_from_electric_induction_divergence as gauss_law,)

from symplyphysics.core.experimental.coordinate_systems import (CoordinateVector, CARTESIAN,
    AppliedPoint)

# The electric induction vector field is given. For its distribution, the electric intensity is
# known, which is equal to 10 kilovolt per meter. A point in space in a three-dimensional
# coordinate system is also given: the x coordinate is 0.493 meter, the y coordinate is 0.5 meter,
# the z coordinate is 0.01249 meter. Then the bulk density of the electric charge at a given point
# will be equal to 244.24 [nanocoulomb / meter^3].

Args = namedtuple("Args", ["electrical_induction", "point"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    electric_intensity = Quantity(10 * units.kilo * units.volt / units.meter)

    x, y, z = CARTESIAN.base_scalars

    electrical_induction = CoordinateVector([
        sin(x / Quantity(1 * units.meter)) * electric_intensity * quantities.vacuum_permittivity,
        sin(y / Quantity(1 * units.meter)) * electric_intensity * quantities.vacuum_permittivity,
        sin(z / Quantity(1 * units.meter)) * electric_intensity * quantities.vacuum_permittivity
    ], CARTESIAN)

    point = AppliedPoint([
        Quantity(0.493 * units.meter),
        Quantity(0.5 * units.meter),
        Quantity(0.01249 * units.meter),
    ], CARTESIAN)

    return Args(electrical_induction=electrical_induction, point=point)


def test_basic_bulk_density_of_electric_charge(test_args: Args) -> None:
    result = gauss_law.calculate_charge_volumetric_density_at_point(test_args.electrical_induction,
        test_args.point)

    assert_equal(result, 244.24 * units.nano * units.coulomb / units.meter**3)


def test_bad_coordinates(test_args: Args) -> None:
    bad_coordinate = Quantity(1 * units.coulomb)

    with raises(errors.UnitsError):
        point = AppliedPoint([bad_coordinate, 0, 0], CARTESIAN)
        gauss_law.calculate_charge_volumetric_density_at_point(test_args.electrical_induction,
            point)
    with raises(TypeError):
        point = AppliedPoint([100, 0, 0], CARTESIAN)
        gauss_law.calculate_charge_volumetric_density_at_point(test_args.electrical_induction,
            point)

    with raises(errors.UnitsError):
        point = AppliedPoint([0, bad_coordinate, 0], CARTESIAN)
        gauss_law.calculate_charge_volumetric_density_at_point(test_args.electrical_induction,
            point)
    with raises(TypeError):
        point = AppliedPoint([0, 100, 0], CARTESIAN)
        gauss_law.calculate_charge_volumetric_density_at_point(test_args.electrical_induction,
            point)

    with raises(errors.UnitsError):
        point = AppliedPoint([0, 0, bad_coordinate], CARTESIAN)
        gauss_law.calculate_charge_volumetric_density_at_point(test_args.electrical_induction,
            point)
    with raises(TypeError):
        point = AppliedPoint([0, 0, 100], CARTESIAN)
        gauss_law.calculate_charge_volumetric_density_at_point(test_args.electrical_induction,
            point)
