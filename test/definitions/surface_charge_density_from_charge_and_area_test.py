from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_approx,
    errors,
    units,
    Quantity,
    SI,
    convert_to,
)

from symplyphysics.definitions import surface_charge_density_from_charge_and_area as surface_charge_density


@fixture(name="test_args")
def test_args_fixture():
    charge = Quantity(20 * units.coulomb)
    area = Quantity(4 * units.meter**2)
    Args = namedtuple("Args", ["charge", "area"])
    return Args(
        charge=charge,
        area=area,
    )


def test_basic_surface_charge_density(test_args):
    result = surface_charge_density.calculate_surface_charge_density(test_args.charge,
        test_args.area)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.charge / units.area)
    result_surface_charge_density = convert_to(result, units.coulomb / units.meter**2).evalf(5)
    assert_approx(result_surface_charge_density, 5)


def test_bad_charge(test_args):
    bad_charge = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        surface_charge_density.calculate_surface_charge_density(bad_charge, test_args.area)
    with raises(TypeError):
        surface_charge_density.calculate_surface_charge_density(100, test_args.area)


def test_bad_area(test_args):
    bad_area = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        surface_charge_density.calculate_surface_charge_density(test_args.charge, bad_area)
    with raises(TypeError):
        surface_charge_density.calculate_surface_charge_density(test_args.charge, 100)
