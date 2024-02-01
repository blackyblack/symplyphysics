from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    errors,
    units,
    Quantity,
    SI,
    convert_to,
)

from symplyphysics.laws.thermodynamics import average_kinetic_energy_of_molecules_from_temperature as average_kinetic_energy


@fixture(name="test_args")
def test_args_fixture():
    temperature = Quantity(300 * units.kelvin)
    Args = namedtuple("Args", ["temperature"])
    return Args(temperature=temperature)


def test_basic_average_kinetic_energy(test_args):
    result = average_kinetic_energy.calculate_average_kinetic_energy(test_args.temperature)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.energy)
    result_energy = convert_to(result, units.joule).evalf(3)
    assert result_energy == approx(6.21e-21, abs=1e-23)


def test_bad_temperature():
    bad_temperature = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        average_kinetic_energy.calculate_average_kinetic_energy(bad_temperature)
    with raises(TypeError):
        average_kinetic_energy.calculate_average_kinetic_energy(100)
