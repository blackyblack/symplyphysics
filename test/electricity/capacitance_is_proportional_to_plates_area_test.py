from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    errors,
    units,
    convert_to,
    Quantity,
    SI,
)
from symplyphysics.laws.electricity import capacitance_is_proportional_to_plates_area as capacitance_law

# Description
## Assert we have a capacitor with 20mm2 plate area and 5micrometer clearance with fielectric permeability of insulator of 5.
## According to online calculator
## (https://matematika-club.ru/ehlektroemkost-kondensatora-kalkulyator-onlajn)
## its capacitance should be  177.08 pF or 0.000000000177083756352408 Farad.


@fixture(name="test_args")
def test_args_fixture():
    permeability = 5.0
    area = Quantity(20 * (units.milli * units.meter)**2)
    clearance = Quantity(5 * units.micro * units.meter)
    Args = namedtuple("Args", ["permeability", "area", "clearance"])
    return Args(permeability=permeability, area=area, clearance=clearance)


def test_basic_capacitance(test_args):
    result = capacitance_law.calculate_capacitance(test_args.permeability, test_args.area, test_args.clearance)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.capacitance)    
    result_capacitance = convert_to(result, units.farad).subs({units.pico: 1, units.farad: 1}).evalf(8)    
    assert result_capacitance == approx(0.000000000177083756352408, 0.000000000001)
    #assert result_capacitance == approx(177.08, 0.000000000001)


def test_bad_area(test_args):
    ab = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        capacitance_law.calculate_capacitance(test_args.permeability, ab, test_args.clearance)
    with raises(TypeError):
        capacitance_law.calculate_capacitance(test_args.permeability, 100, test_args.clearance)


def test_bad_clearance(test_args):
    cb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        capacitance_law.calculate_capacitance(test_args.permeability, test_args.area, cb)
    with raises(TypeError):
        capacitance_law.calculate_capacitance(test_args.permeability, test_args.area, 100)
