from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, units, Quantity, errors)
from symplyphysics.laws.electricity.chemistry import energy_of_an_electron_in_a_hydrogen_atom as energy_law

# Description
## The energy value is known for the radius of the first Bohr orbit of the electron. It is equal to 13.6 [eV].
## https://en.wikipedia.org/wiki/Bohr_model

Args = namedtuple("Args", ["radius_of_electron"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    radius_of_electron = Quantity(5.291e-11 * units.meter)
    return Args(radius_of_electron=radius_of_electron)


def test_basic_energy_of_electron(test_args: Args) -> None:
    result = energy_law.calculate_energy_of_electron(test_args.radius_of_electron)
    assert_equal(result, 13.6 * units.electronvolt)


def test_bad_radius_of_electron() -> None:
    radius_of_electron = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        energy_law.calculate_energy_of_electron(radius_of_electron)
    with raises(TypeError):
        energy_law.calculate_energy_of_electron(100)
