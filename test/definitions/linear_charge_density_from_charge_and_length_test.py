from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    errors,
    units,
    Quantity,
    SI,
    convert_to,
    dimensionless
)

from symplyphysics.definitions import linear_charge_density_from_charge_and_length as linear_charge_density


@fixture(name="test_args")
def test_args_fixture():
    charge = Quantity(8 * units.coulomb)
    length = Quantity(4 * units.meter)
    Args = namedtuple("Args", ["charge", "length"])
    return Args(
        charge=charge,
        length=length,
               )


def test_basic_linear_charge_density(test_args):
    result = linear_charge_density.calculate_linear_charge_density(test_args.charge, test_args.length)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.charge / units.length)
    result_linear_charge_density = convert_to(result, units.coulomb / units.meter).evalf(5)
    assert result_linear_charge_density == approx(2, 0.001)


def test_bad_charge(test_args):
    bad_charge = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        linear_charge_density.calculate_linear_charge_density(bad_charge, test_args.length)
    with raises(TypeError):
        linear_charge_density.calculate_linear_charge_density(100, test_args.length)


def test_bad_length(test_args):
    bad_length = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        linear_charge_density.calculate_linear_charge_density(test_args.charge, bad_length)
    with raises(TypeError):
        linear_charge_density.calculate_linear_charge_density(test_args.charge, 100)
