from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_approx, errors, units, Quantity, SI, convert_to, prefixes)
from symplyphysics.laws.optics import radiation_intensity_from_energy_area_time as intensity_law

# Description
## Let the energy of the incident light be 110 microjoule, the surface area is 5 [centimeter^2], and the time
## is 5 second. Then the radiation power is 44 [milliwatt / meter^2].
## https://www.indigomath.ru//raschety/lXPL2K.html


@fixture(name="test_args")
def test_args_fixture():
    energy = Quantity(110 * prefixes.micro * units.joule)
    area = Quantity(5 * units.centimeter**2)
    time = Quantity(5 * units.second)
    Args = namedtuple("Args", ["energy", "area", "time"])
    return Args(energy=energy, area=area, time=time)


def test_basic_intensity(test_args):
    result = intensity_law.calculate_intensity(test_args.energy, test_args.area, test_args.time)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.power / units.area)
    result = convert_to(result, prefixes.milli * units.watt / units.meter**2).evalf(2)
    assert_approx(result, 44)


def test_bad_energy(test_args):
    energy = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        intensity_law.calculate_intensity(energy, test_args.area, test_args.time)
    with raises(TypeError):
        intensity_law.calculate_intensity(100, test_args.area, test_args.time)


def test_bad_area(test_args):
    area = Quantity(1 * units.kelvin)
    with raises(errors.UnitsError):
        intensity_law.calculate_intensity(test_args.energy, area, test_args.time)
    with raises(TypeError):
        intensity_law.calculate_intensity(test_args.energy, 100, test_args.time)


def test_bad_time(test_args):
    time = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        intensity_law.calculate_intensity(test_args.energy, test_args.area, time)
    with raises(TypeError):
        intensity_law.calculate_intensity(test_args.energy, test_args.area, 100)
