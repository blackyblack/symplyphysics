from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    units,
    SI,
    convert_to,
    Quantity,
    errors
)
from symplyphysics.laws.condensed_matter import current_density_from_concentration_and_velocity_of_charge_carriers as current_law

# Description
## Calculate the current density for a copper conductor. It is known that it is equal to 1e6 [amper / meter**2]
## https://online.mephi.ru/courses/physics/electricity/data/course/4/4.1.html


@fixture(name="test_args")
def test_args_fixture():
    charge_carriers_concentration = Quantity(8.4e28 * (1 / units.meter**3))
    drift_velocity = Quantity(7.4e-5 * (units.meter / units.second))
    charge_electron = Quantity(1.6e-19 * units.coulomb)

    Args = namedtuple("Args", ["charge_carriers_concentration", "drift_velocity", "charge_electron"])
    return Args(
        charge_carriers_concentration=charge_carriers_concentration,
        drift_velocity=drift_velocity, charge_electron=charge_electron)


def test_basic_current(test_args):
    result = current_law.calculate_current(test_args.charge_carriers_concentration, test_args.drift_velocity, test_args.charge_electron)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.current / units.area)
    result = convert_to(result, units.ampere / units.meter**2).evalf(5)
    assert result == approx(1e6, rel=0.01)


def test_bad_charge_carriers_concentration(test_args):
    charge_carriers_concentration = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        current_law.calculate_current(charge_carriers_concentration, test_args.drift_velocity, test_args.charge_electron)
    with raises(TypeError):
        current_law.calculate_current(100, test_args.drift_velocity, test_args.charge_electron)


def test_bad_drift_velocity(test_args):
    drift_velocity = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        current_law.calculate_current(test_args.charge_carriers_concentration, drift_velocity, test_args.charge_electron)
    with raises(TypeError):
        current_law.calculate_current(test_args.charge_carriers_concentration, 100, test_args.charge_electron)


def test_bad_charge_electron(test_args):
    charge_electron = Quantity(1 * units.kelvin)
    with raises(errors.UnitsError):
        current_law.calculate_current(test_args.charge_carriers_concentration, test_args.drift_velocity, charge_electron)
    with raises(TypeError):
        current_law.calculate_current(test_args.charge_carriers_concentration, test_args.drift_velocity, 100)
