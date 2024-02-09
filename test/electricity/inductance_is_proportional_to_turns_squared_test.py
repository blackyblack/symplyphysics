from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, errors, units, Quantity, prefixes)
from symplyphysics.laws.electricity import inductance_is_proportional_to_turns_squared as coil_law

# Description
## Assert we have a capacitor with 100 turns on a core with relative permeability of 2, length of 5 mm and turn area is 0.002 m^2.
## According to online calculator (https://goo.su/WBwDZXS)
## its inductance should be 10.05309649 mH.

Args = namedtuple("Args", ["permeability", "turncount", "turn_area", "length"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    permeability = 2
    turncount = 100
    turn_area = Quantity(0.002 * (units.meter)**2)
    length = Quantity(5 * prefixes.milli * units.meter)
    return Args(permeability=permeability, turncount=turncount, turn_area=turn_area, length=length)


def test_basic_inductance(test_args: Args) -> None:
    result = coil_law.calculate_inductance(test_args.permeability, test_args.turncount,
        test_args.turn_area, test_args.length)
    assert_equal(result, 10.053 * prefixes.milli * units.henry)


def test_bad_area(test_args: Args) -> None:
    ab = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        coil_law.calculate_inductance(test_args.permeability, test_args.turncount, ab,
            test_args.length)
    with raises(TypeError):
        coil_law.calculate_inductance(test_args.permeability, test_args.turncount, 100,
            test_args.length)


def test_bad_length(test_args: Args) -> None:
    lb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        coil_law.calculate_inductance(test_args.permeability, test_args.turncount,
            test_args.turn_area, lb)
    with raises(TypeError):
        coil_law.calculate_inductance(test_args.permeability, test_args.turncount,
            test_args.turn_area, 100)
