from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    units,
    errors,
    Quantity,
)
from symplyphysics.laws.condensed_matter import equilibrium_voltage_difference_in_pn_junction_via_concentrations as barrier_law

# Description
## For concentration of donors equal to 3e16 [1 / cm^3] and concentration of acceptors equal to 2e12 [1 / cm^3]
## the height of the potential barrier is known. It is equal to 0.529 volts. You can check it on the website below:
## https://cleanroom.byu.edu/pn_junction
## It is known that the concentration of intrinsic charge carriers in silicon is equal to 1e10 [1 / cm^3].
## https://www.ioffe.ru/SVA/NSM/Semicond/Si/bandstr.html

Args = namedtuple("Args", [
    "donors_concentration", "acceptors_concentration", "charge_carriers_concentration",
    "temperature", "charge_electron"
])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    donors_concentration = Quantity(3e16 * (1 / units.centimeter**3))
    acceptors_concentration = Quantity(2e12 * (1 / units.centimeter**3))
    charge_carriers_concentration = Quantity(1e10 * (1 / units.centimeter**3))
    temperature = Quantity(300 * units.kelvin)
    charge_electron = Quantity(1.6e-19 * units.coulomb)

    return Args(donors_concentration=donors_concentration,
        acceptors_concentration=acceptors_concentration,
        charge_carriers_concentration=charge_carriers_concentration,
        temperature=temperature,
        charge_electron=charge_electron)


def test_basic_height_barrier(test_args: Args) -> None:
    result = barrier_law.calculate_height_barrier(test_args.donors_concentration,
        test_args.acceptors_concentration, test_args.charge_carriers_concentration,
        test_args.temperature, test_args.charge_electron)
    # NOTE: onsite intrinsic charge carriers are caclulated instead of hardcoded 1e10 value,
    #       therefore result in our test is not very accurate
    assert_equal(result, 0.529 * units.volt, relative_tolerance=0.1)


def test_bad_donors_concentration(test_args: Args) -> None:
    donors_concentration = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        barrier_law.calculate_height_barrier(donors_concentration,
            test_args.acceptors_concentration, test_args.charge_carriers_concentration,
            test_args.temperature, test_args.charge_electron)
    with raises(TypeError):
        barrier_law.calculate_height_barrier(100, test_args.acceptors_concentration,
            test_args.charge_carriers_concentration, test_args.temperature,
            test_args.charge_electron)


def test_bad_acceptors_concentration(test_args: Args) -> None:
    acceptors_concentration = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        barrier_law.calculate_height_barrier(test_args.donors_concentration,
            acceptors_concentration, test_args.charge_carriers_concentration, test_args.temperature,
            test_args.charge_electron)
    with raises(TypeError):
        barrier_law.calculate_height_barrier(test_args.donors_concentration, 100,
            test_args.charge_carriers_concentration, test_args.temperature,
            test_args.charge_electron)


def test_bad_charge_carriers_concentration(test_args: Args) -> None:
    charge_carriers_concentration = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        barrier_law.calculate_height_barrier(test_args.donors_concentration,
            test_args.acceptors_concentration, charge_carriers_concentration, test_args.temperature,
            test_args.charge_electron)
    with raises(TypeError):
        barrier_law.calculate_height_barrier(test_args.donors_concentration,
            test_args.acceptors_concentration, 100, test_args.temperature,
            test_args.charge_electron)


def test_bad_temperature(test_args: Args) -> None:
    temperature = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        barrier_law.calculate_height_barrier(test_args.donors_concentration,
            test_args.acceptors_concentration, test_args.charge_carriers_concentration, temperature,
            test_args.charge_electron)
    with raises(TypeError):
        barrier_law.calculate_height_barrier(test_args.donors_concentration,
            test_args.acceptors_concentration, test_args.charge_carriers_concentration, 100,
            test_args.charge_electron)


def test_bad_charge_electron(test_args: Args) -> None:
    charge_electron = Quantity(1 * units.kelvin)
    with raises(errors.UnitsError):
        barrier_law.calculate_height_barrier(test_args.donors_concentration,
            test_args.acceptors_concentration, test_args.charge_carriers_concentration,
            test_args.temperature, charge_electron)
    with raises(TypeError):
        barrier_law.calculate_height_barrier(test_args.donors_concentration,
            test_args.acceptors_concentration, test_args.charge_carriers_concentration,
            test_args.temperature, 100)
