from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.hydro.hydrostatic_pressure_via_density_and_height import calculate_hydrostatic_pressure

# Description
# A body at a depth of 10 m in water (1000 kg/m³) with a free fall acceleration of 9.81 m/s²
# should experience a hydrostatic pressure of 98100 Pa.

Args = namedtuple("Args", ["rho", "h"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    rho = Quantity(1000 * units.kilogram / units.meter**3)
    h = Quantity(10 * units.meter)
    return Args(rho=rho, h=h)


def test_hydrostatic_pressure(test_args: Args) -> None:
    result = calculate_hydrostatic_pressure(test_args.rho, test_args.h)
    assert_equal(result, 98100 * units.pascal)


def test_bad_density(test_args: Args) -> None:
    bad_density = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        calculate_hydrostatic_pressure(bad_density, test_args.h)
    with raises(TypeError):
        calculate_hydrostatic_pressure(100, test_args.h)


def test_bad_depth(test_args: Args) -> None:
    bad_depth = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        calculate_hydrostatic_pressure(test_args.rho, bad_depth)
    with raises(TypeError):
        calculate_hydrostatic_pressure(test_args.rho, 100)
