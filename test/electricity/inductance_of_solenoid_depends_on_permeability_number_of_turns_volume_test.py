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
from symplyphysics.laws.electricity import inductance_of_solenoid_depends_on_permeability_number_of_turns_volume as inductance_law

# Description
## Number of turns per unit length is 10, magnetic permeability of solenoid core is 10,
## and volume of solenoid is 0.1 [meter^3]. Then the inductance of the solenoid is 125.6 microhenry.
## https://physics.icalculator.com/self-inductance-calculator.html


@fixture(name="test_args")
def test_args_fixture():
    relative_permeability = 10
    number_of_turns = Quantity(10 * (1 / units.meter))
    volume = Quantity(0.1 * units.meter**3)

    Args = namedtuple("Args", ["relative_permeability", "number_of_turns", "volume"])
    return Args(relative_permeability=relative_permeability, number_of_turns=number_of_turns, volume=volume)


def test_basic_inductance(test_args):
    result = inductance_law.calculate_inductance(test_args.relative_permeability, test_args.number_of_turns, test_args.volume)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.inductance)
    result_voltage = convert_to(result, prefixes.micro * units.henry).evalf(5)
    assert result_voltage == approx(125.6, 0.01)


def test_bad_relative_permeability(test_args):
    relative_permeability = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        inductance_law.calculate_inductance(relative_permeability, test_args.number_of_turns, test_args.volume)


def test_bad_number_of_turns(test_args):
    number_of_turns = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        inductance_law.calculate_inductance(test_args.relative_permeability, number_of_turns, test_args.volume)
    with raises(TypeError):
        inductance_law.calculate_inductance(test_args.relative_permeability, 100, test_args.volume)


def test_bad_volume(test_args):
    volume = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        inductance_law.calculate_inductance(test_args.relative_permeability, test_args.number_of_turns, volume)
    with raises(TypeError):
        inductance_law.calculate_inductance(test_args.relative_permeability, test_args.number_of_turns, 100)
