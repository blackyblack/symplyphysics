from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (errors, units, Quantity, SI, convert_to)

from symplyphysics.laws.electricity import volume_charge_density_from_charge_and_volume as volume_charge_density


@fixture(name="test_args")
def test_args_fixture():
    charge = Quantity(6 * units.coulomb)
    volume = Quantity(3 * units.meter**3)
    Args = namedtuple("Args", ["charge", "volume"])
    return Args(
        charge=charge,
        volume=volume,
    )


def test_basic_volume_charge_density(test_args):
    result = volume_charge_density.calculate_volume_charge_density(test_args.charge,
        test_args.volume)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.charge / units.volume)
    result_volume_charge_density = convert_to(result, units.coulomb / units.meter**3).evalf(5)
    assert result_volume_charge_density == approx(2, 0.001)


def test_bad_charge(test_args):
    bad_charge = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        volume_charge_density.calculate_volume_charge_density(bad_charge, test_args.volume)
    with raises(TypeError):
        volume_charge_density.calculate_volume_charge_density(100, test_args.volume)


def test_bad_volume(test_args):
    bad_volume = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        volume_charge_density.calculate_volume_charge_density(test_args.charge, bad_volume)
    with raises(TypeError):
        volume_charge_density.calculate_volume_charge_density(test_args.charge, 100)
