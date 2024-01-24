from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (units, SI, convert_to, Quantity, errors)
from symplyphysics.laws.condensed_matter import drift_velocity_of_charge_carriers as velocity_law

# Description
## At high values of the field, the drift velocity of electrons in silicon is approximately 1e5 [meter / second].
## https://foez.narod.ru/19.htm#:~:text=%D0%92%20%D0%BA%D1%80%D0%B5%D0%BC%D0%BD%D0%B8%D0%B8%20%D0%BF%D1%80%D0%B8%20%D0%A2%3D300,%C3%97104%20%D0%BC%2Fc.


@fixture(name="test_args")
def test_args_fixture():
    charge_carriers_mobility = Quantity(0.17 * (units.meter**2 / units.volt / units.second))
    electric_intensity = Quantity(5.5e5 * (units.volt / units.meter))

    Args = namedtuple("Args", ["charge_carriers_mobility", "electric_intensity"])
    return Args(charge_carriers_mobility=charge_carriers_mobility,
        electric_intensity=electric_intensity)


def test_basic_drift_velocity(test_args):
    result = velocity_law.calculate_velocity(test_args.charge_carriers_mobility,
        test_args.electric_intensity)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.length / units.time)
    result = convert_to(result, units.meter / units.second).evalf(5)
    assert result == approx(1e5, rel=0.1)


def test_bad_charge_carriers_mobility(test_args):
    charge_carriers_mobility = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        velocity_law.calculate_velocity(charge_carriers_mobility, test_args.electric_intensity)
    with raises(TypeError):
        velocity_law.calculate_velocity(100, test_args.electric_intensity)


def test_bad_electric_intensity(test_args):
    electric_intensity = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        velocity_law.calculate_velocity(test_args.charge_carriers_mobility, electric_intensity)
    with raises(TypeError):
        velocity_law.calculate_velocity(test_args.charge_carriers_mobility, 100)
