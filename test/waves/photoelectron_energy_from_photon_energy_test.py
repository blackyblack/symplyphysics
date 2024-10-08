from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.waves import photoelectron_energy_from_photon_energy as photoeffect

# Description
## The metal surface with the work function 3.81e-19 J is irradiated with photons of frequency 9.23e14 Hz.
## The maximum kinetic energy of the photoelectrons will be 2.3e-19 J
## Source: https://zaochnik.ru/blog/zadachi-na-temu-fotony-i-fotoeffekt-s-resheniem/?ysclid=lqmlagmvdp85340848 (4)

Args = namedtuple("Args", ["energy", "work"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    frequency = Quantity(9.224e14 * units.hertz)
    energy = units.planck * frequency
    work = Quantity(3.81e-19 * units.joule)
    return Args(energy=energy, work=work)


def test_basic_energy(test_args: Args) -> None:
    result = photoeffect.calculate_max_kinetic_energy(test_args.energy, test_args.work)
    assert_equal(result, 2.3e-19 * units.joule)


def test_bad_energy(test_args: Args) -> None:
    eb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        photoeffect.calculate_max_kinetic_energy(eb, test_args.work)
    with raises(TypeError):
        photoeffect.calculate_max_kinetic_energy(100, test_args.work)


def test_bad_work_function(test_args: Args) -> None:
    wb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        photoeffect.calculate_max_kinetic_energy(test_args.energy, wb)
    with raises(TypeError):
        photoeffect.calculate_max_kinetic_energy(test_args.energy, 100)
