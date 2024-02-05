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

from symplyphysics.laws.thermodynamics import pressure_from_the_concentration_of_molecules_and_kinetic_energy as ideal_gas_pressure


@fixture(name="test_args")
def test_args_fixture():
    average_kinetic_energy = Quantity(2.4e-20 * units.joule)
    molecules_concentration = Quantity(5e24 * (1 / units.meter**3))
    Args = namedtuple("Args", ["average_kinetic_energy", "molecules_concentration"])
    return Args(average_kinetic_energy=average_kinetic_energy,
        molecules_concentration=molecules_concentration)


def test_basic_pressure(test_args):
    result = ideal_gas_pressure.calculate_pressure(test_args.molecules_concentration,
        test_args.average_kinetic_energy)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.pressure)
    result_pressure = convert_to(result, units.pascal).evalf(3)
    assert_approx(result_pressure, 80000)


def test_bad_energy(test_args):
    bad_energy = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        ideal_gas_pressure.calculate_pressure(test_args.molecules_concentration, bad_energy)
    with raises(TypeError):
        ideal_gas_pressure.calculate_pressure(test_args.molecules_concentration, 100)


def test_bad_concentration(test_args):
    bad_concentration = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        ideal_gas_pressure.calculate_pressure(bad_concentration, test_args.average_kinetic_energy)
    with raises(TypeError):
        ideal_gas_pressure.calculate_pressure(100, test_args.average_kinetic_energy)
