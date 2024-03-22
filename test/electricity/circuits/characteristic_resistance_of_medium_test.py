from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, units, Quantity, errors)
from symplyphysics.laws.electricity.circuits import characteristic_resistance_of_medium as resistance_law

# Description
## The magnetic permeability of the medium is 1, the relative_permittivity is 2.2. Then the characteristic
## resistance is 254.167 ohm.

Args = namedtuple("Args", ["relative_permittivity", "relative_permeability"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    relative_permittivity = 2.2
    relative_permeability = 1

    return Args(relative_permittivity=relative_permittivity, relative_permeability=relative_permeability)


def test_basic_resistance(test_args: Args) -> None:
    result = resistance_law.calculate_resistance(test_args.relative_permittivity, test_args.relative_permeability)
    assert_equal(result, 254.167 * units.ohm)


def test_bad_relative_permittivity(test_args: Args) -> None:
    relative_permittivity = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        resistance_law.calculate_resistance(relative_permittivity, test_args.relative_permeability)


def test_bad_relative_permeability(test_args: Args) -> None:
    relative_permeability = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        resistance_law.calculate_resistance(test_args.relative_permittivity, relative_permeability)
