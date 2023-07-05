from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (errors, units, Quantity, SI, convert_to)
from symplyphysics.laws.electricity import inductance_is_proportional_to_turns_squared as coil_law

# Description
## Assert we have a capacitor with 100 turns on a core with relative permeability of 2, length of 5mm and turn area is 0.002m2.
## According to online calculator (https://www.asutpp.ru/onlayn-kalkulyator-rascheta-induktivnosti-katushki-s-vozdushnym-serdechnikom.html)
## its inductance should be  10.555751316062 mH.


@fixture(name="test_args")
def test_args_fixture():
    permeability = 2
    turncount = 100   
    turn_area = Quantity(0.002 * (units.meter)**2)
    length = Quantity(5 * units.milli * units.meter)
    Args = namedtuple("Args", ["permeability", "turncount", "turn_area", "length"])
    return Args(permeability=permeability, turncount=turncount, turn_area=turn_area, length=length)


def test_basic_inductance(test_args):    
    result = coil_law.calculate_inductance(test_args.permeability, test_args.turncount, test_args.turn_area, test_args.length)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.inductance)
    result_inductance = convert_to(result, units.milli * units.henry).evalf(6)
    assert result_inductance == approx(10.555751316062, 0.0001)

'''
def test_bad_area(test_args):
    ab = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        capacitance_law.calculate_capacitance(test_args.permeability, ab, test_args.distance)
    with raises(TypeError):
        capacitance_law.calculate_capacitance(test_args.permeability, 100, test_args.distance)


def test_bad_distance(test_args):
    db = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        capacitance_law.calculate_capacitance(test_args.permeability, test_args.area, db)
    with raises(TypeError):
        capacitance_law.calculate_capacitance(test_args.permeability, test_args.area, 100)
'''
