from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    errors,
    units,
    convert_to,
    Quantity,
    SI,
)
from symplyphysics.laws.kinematic import maximum_height_from_velocity as maximum_height_law

# Description
## If a body is thrown vertically upwards at a speed of 10 m/s, it will rise to a height of  (10*10) /2*9,8 = 5.1 meters

@fixture(name="test_args")
def test_args_fixture():
    v = Quantity(10 * units.meter / units.second)
    Args = namedtuple("Args", ["v"])
    return Args(v=v)

def test_basic_height(test_args):
    result = maximum_height_law.calculate_height(test_args.v)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.length)
    result_vector = convert_to(result, units.meter).evalf(2)
    assert result_vector == approx(5.10, 0.01)

def test_bad_velocity(test_args):
    vb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
       maximum_height_law.calculate_height(vb)
    with raises(TypeError):
       maximum_height_law.calculate_height(100)