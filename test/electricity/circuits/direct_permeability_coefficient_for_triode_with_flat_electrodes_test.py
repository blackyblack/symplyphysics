from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, units, Quantity, errors)
from symplyphysics.laws.electricity.circuits import direct_permeability_coefficient_for_triode_with_flat_electrodes as coefficient_law

# Description
## The grid step is equal to 3 millimeter. The distance from the anode to the grid is 5 centimeter.
## The first tabular coefficient is 2. The second tabular coefficient is 0.1.
## Then the direct permeability coefficient of the grid for a triode with flat electrodes will be 1.2.

Args = namedtuple("Args", ["tabular_coefficient_T", "grid_step", "distance_to_grid", "tabular_coefficient_d"])

@fixture(name="test_args")
def test_args_fixture() -> Args:
    tabular_coefficient_T = 2
    grid_step = Quantity(3 * units.millimeter)
    distance_to_grid = Quantity(5 * units.centimeter)
    tabular_coefficient_d = 0.1

    return Args(tabular_coefficient_T=tabular_coefficient_T,
                grid_step=grid_step,
                distance_to_grid=distance_to_grid,
                tabular_coefficient_d=tabular_coefficient_d)


def test_basic_direct_permeability_coefficient(test_args: Args) -> None:
    result = coefficient_law.calculate_direct_permeability_coefficient(
        test_args.tabular_coefficient_T, test_args.grid_step, test_args.distance_to_grid, test_args.tabular_coefficient_d)
    assert_equal(result, 1.2)


def test_bad_tabular_coefficient_T(test_args: Args) -> None:
    tabular_coefficient_T = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        coefficient_law.calculate_direct_permeability_coefficient(tabular_coefficient_T, test_args.grid_step, test_args.distance_to_grid, test_args.tabular_coefficient_d)


def test_bad_grid_step(test_args: Args) -> None:
    grid_step = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        coefficient_law.calculate_direct_permeability_coefficient(test_args.tabular_coefficient_T, grid_step, test_args.distance_to_grid, test_args.tabular_coefficient_d)
    with raises(TypeError):
        coefficient_law.calculate_direct_permeability_coefficient(test_args.tabular_coefficient_T, 100, test_args.distance_to_grid, test_args.tabular_coefficient_d)


def test_bad_distance_to_grid(test_args: Args) -> None:
    distance_to_grid = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        coefficient_law.calculate_direct_permeability_coefficient(test_args.tabular_coefficient_T, test_args.grid_step, distance_to_grid, test_args.tabular_coefficient_d)
    with raises(TypeError):
        coefficient_law.calculate_direct_permeability_coefficient(test_args.tabular_coefficient_T, test_args.grid_step, 100, test_args.tabular_coefficient_d)


def test_bad_tabular_coefficient_d(test_args: Args) -> None:
    tabular_coefficient_d = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        coefficient_law.calculate_direct_permeability_coefficient(test_args.tabular_coefficient_T, test_args.grid_step, test_args.distance_to_grid, tabular_coefficient_d)
