from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_approx,
    errors,
    units,
    convert_to,
    Quantity,
    SI,
)
from symplyphysics.laws.waves import photoelectron_energy_from_frequency as photoeffect

# Description
## The metal surface with the work function 3.81e-19 J is irradiated with photons of frequency 9.23e14 Hz.
## The maximum kinetic energy of the photoelectrons will be 2.3e-19 J
## Source: https://zaochnik.ru/blog/zadachi-na-temu-fotony-i-fotoeffekt-s-resheniem/?ysclid=lqmlagmvdp85340848 (4)


@fixture(name="test_args")
def test_args_fixture():
    frequency = Quantity(9.224e14 * units.hertz)
    work = Quantity(3.81e-19 * units.joule)
    Args = namedtuple("Args", ["frequency", "work"])
    return Args(frequency=frequency, work=work)


def test_basic_energy(test_args):
    result = photoeffect.calculate_max_kinetic_energy(test_args.frequency, test_args.work)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.energy)
    result_energy = convert_to(result, units.joule).evalf(4)
    assert_approx(result_energy, 2.3e-19)


def test_bad_frequency(test_args):
    fb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        photoeffect.calculate_max_kinetic_energy(fb, test_args.work)
    with raises(TypeError):
        photoeffect.calculate_max_kinetic_energy(100, test_args.work)


def test_bad_work_function(test_args):
    wb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        photoeffect.calculate_max_kinetic_energy(test_args.frequency, wb)
    with raises(TypeError):
        photoeffect.calculate_max_kinetic_energy(test_args.frequency, 100)
