from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
    prefixes,
)
from symplyphysics.laws.electricity import inductance_of_solenoid_depends_on_permeability_number_of_turns_volume as inductance_law

# Description
## Number of turns per unit length is 10, magnetic permeability of solenoid core is 10,
## and volume of solenoid is 0.1 [meter^3]. Then the inductance of the solenoid is 125.6 microhenry.
## https://physics.icalculator.com/self-inductance-calculator.html

Args = namedtuple("Args", ["relative_permeability", "number_of_turns", "volume"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    relative_permeability = 10
    number_of_turns = Quantity(10 * (1 / units.meter))
    volume = Quantity(0.1 * units.meter**3)
    return Args(relative_permeability=relative_permeability,
        number_of_turns=number_of_turns,
        volume=volume)


def test_basic_inductance(test_args: Args) -> None:
    result = inductance_law.calculate_inductance(test_args.relative_permeability,
        test_args.number_of_turns, test_args.volume)
    assert_equal(result, 125.6 * prefixes.micro * units.henry)


def test_bad_relative_permeability(test_args: Args) -> None:
    relative_permeability = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        inductance_law.calculate_inductance(relative_permeability, test_args.number_of_turns,
            test_args.volume)


def test_bad_number_of_turns(test_args: Args) -> None:
    number_of_turns = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        inductance_law.calculate_inductance(test_args.relative_permeability, number_of_turns,
            test_args.volume)
    with raises(TypeError):
        inductance_law.calculate_inductance(test_args.relative_permeability, 100, test_args.volume)


def test_bad_volume(test_args: Args) -> None:
    volume = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        inductance_law.calculate_inductance(test_args.relative_permeability,
            test_args.number_of_turns, volume)
    with raises(TypeError):
        inductance_law.calculate_inductance(test_args.relative_permeability,
            test_args.number_of_turns, 100)
