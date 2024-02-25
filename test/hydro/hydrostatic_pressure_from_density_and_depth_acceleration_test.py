from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.hydro.hydrostatic_pressure_from_density_and_depth_acceleration import calculate_hydrostatic_pressure

# Description
# A body at a depth of 10 m in water (1000 kg/m³) with a free fall acceleration of 9.81 m/s²
# should experience a hydrostatic pressure of 98100 Pa.

Args = namedtuple("Args", ["rho", "h", "acceleration"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    rho = Quantity(1000 * units.kilogram / units.meter**3)
    h = Quantity(10 * units.meter)
    acceleration = Quantity(9.81 * (units.meter / units.second**2))
    return Args(rho=rho, h=h, acceleration=acceleration)


def test_hydrostatic_pressure(test_args: Args) -> None:
    result = calculate_hydrostatic_pressure(test_args.rho, test_args.h, test_args.acceleration)
    assert_equal(result, 98100 * units.pascal)


def test_bad_density(test_args: Args) -> None:
    bad_density = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        calculate_hydrostatic_pressure(bad_density, test_args.h, test_args.acceleration)
    with raises(errors.UnitsError):
        calculate_hydrostatic_pressure("not a quantity", test_args.h, test_args.acceleration)


def test_bad_depth(test_args: Args) -> None:
    bad_depth = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        calculate_hydrostatic_pressure(test_args.rho, bad_depth, test_args.acceleration)
    with raises(errors.UnitsError):
        calculate_hydrostatic_pressure(test_args.rho, "not a quantity", test_args.acceleration)


def test_bad_acceleration(test_args: Args) -> None:
    bad_acceleration = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        calculate_hydrostatic_pressure(test_args.rho, test_args.h, bad_acceleration)
    with raises(errors.UnitsError):
        calculate_hydrostatic_pressure(test_args.rho, test_args.h, "not a quantity")
