from collections import namedtuple
from pytest import approx, fixture, raises
from sympy.physics.units import prefixes
from symplyphysics import (
    errors,
    units,
    convert_to,
    Quantity,
    SI,
)
from symplyphysics.laws.electricity import energy_of_magnetic_field_of_solenoid as energy_law

# Description
## Magnetic permeability of solenoid core is 120, intensity of magnetic field is 5 ampere per meter,
## the volume of the solenoid is 0.1 [meter^3]. Then energy of solenoid will be equal to 190 microjoule.
## https://www.indigomath.ru//raschety/Pwhwei.html


@fixture(name="test_args")
def test_args_fixture():
    relative_permeability = 120
    intensity = Quantity(5 * (units.ampere / units.meter))
    volume = Quantity(0.1 * units.meter**3)

    Args = namedtuple("Args", ["relative_permeability", "intensity", "volume"])
    return Args(relative_permeability=relative_permeability, intensity=intensity, volume=volume)


def test_basic_energy(test_args):
    result = energy_law.calculate_energy(test_args.relative_permeability, test_args.intensity, test_args.volume)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.energy)
    result_voltage = convert_to(result, prefixes.micro * units.joule).evalf(5)
    assert result_voltage == approx(190, 0.01)


def test_bad_relative_permeability(test_args):
    relative_permeability = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        energy_law.calculate_energy(relative_permeability, test_args.intensity, test_args.volume)


def test_bad_intensity(test_args):
    intensity = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        energy_law.calculate_energy(test_args.relative_permeability, intensity, test_args.volume)
    with raises(TypeError):
        energy_law.calculate_energy(test_args.relative_permeability, 100, test_args.volume)


def test_bad_volume(test_args):
    volume = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        energy_law.calculate_energy(test_args.relative_permeability, test_args.intensity, volume)
    with raises(TypeError):
        energy_law.calculate_energy(test_args.relative_permeability, test_args.intensity, 100)
