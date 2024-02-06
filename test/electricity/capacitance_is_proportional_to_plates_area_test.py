from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_approx, errors, units, Quantity, SI, convert_to, prefixes)
from symplyphysics.laws.electricity import capacitance_is_proportional_to_plates_area as capacitance_law

# Description
## Assert we have a capacitor with 20mm2 plate area and 5micrometer clearance with dielectric permeability of insulator of 5.
## According to online calculator
## (https://matematika-club.ru/ehlektroemkost-kondensatora-kalkulyator-onlajn)
## its capacitance should be  177.08 pF or 0.000000000177083756352408 Farad.


@fixture(name="test_args")
def test_args_fixture():
    permeability = 5.0
    area = Quantity(20 * (prefixes.milli * units.meter)**2)
    distance = Quantity(5 * prefixes.micro * units.meter)
    Args = namedtuple("Args", ["permeability", "area", "distance"])
    return Args(permeability=permeability, area=area, distance=distance)


def test_basic_capacitance(test_args):
    result = capacitance_law.calculate_capacitance(test_args.permeability, test_args.area,
        test_args.distance)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.capacitance)
    result_capacitance = convert_to(result, prefixes.pico * units.farad).evalf(4)
    assert_approx(result_capacitance, 177.08)


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
