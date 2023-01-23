from collections import namedtuple
from pytest import approx, fixture, raises

from symplyphysics import (
    units, convert_to, SI, errors
)
from symplyphysics.laws.hydro import velocity_from_height as torrichellis_formula

# Description
## In the bottom of water tank there is a small hole. Heigth of liquid in tank is 3m. Velocity of jet according to TOrrichelli's formula should be
## 7.67m/s

@fixture
def test_args():
    g = units.Quantity('g')
    SI.set_quantity_dimension(g, units.length / units.time**2)
    SI.set_quantity_scale_factor(g, 9.80665 * units.meter / units.second**2)

    h = units.Quantity('h')
    SI.set_quantity_dimension(h, units.length)
    SI.set_quantity_scale_factor(h, 3 * units.meter)

    Args = namedtuple('Args', ['g', 'h'])
    return Args(g = g, h = h)

def test_basic_velocity(test_args):
    result = torrichellis_formula.calculate_velocity(test_args.g, test_args.h)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.velocity)
    result_velocity = convert_to(result, units.meter / units.second).subs({units.meter: 1, units.second: 1}).evalf(5)
    assert result_velocity == approx(7.67, 0.0001)

def test_bad_gravity(test_args):
    gb = units.Quantity('gb')
    SI.set_quantity_dimension(gb, units.length)
    SI.set_quantity_scale_factor(gb, 1 * units.meter)

    with raises(errors.UnitsError):
        torrichellis_formula.calculate_velocity(gb, test_args.h)

    with raises(TypeError):
        torrichellis_formula.calculate_velocity(100, test_args.h)

def test_bad_height(test_args):
    hb = units.Quantity('hb')
    SI.set_quantity_dimension(hb, units.charge)
    SI.set_quantity_scale_factor(hb, 1 * units.coulomb)

    with raises(errors.UnitsError):
        torrichellis_formula.calculate_velocity(test_args.g, hb)

    with raises(TypeError):
        torrichellis_formula.calculate_velocity(test_args.g, 100)
        