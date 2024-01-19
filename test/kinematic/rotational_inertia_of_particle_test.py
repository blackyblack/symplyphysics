from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    errors,
    units,
    convert_to,
    Quantity,
    SI,
)
from symplyphysics.laws.kinematic import rotational_inertia_of_particle as rotational_inertia_def

# Description
## A particle of mass m = 1 g rotates around an axis at a radius r = 3 m. Its rotational inertia
## should amount to 0.009 kg*m^2.


@fixture(name="test_args")
def test_args_fixture():
    m = Quantity(1.0 * units.gram)
    r = Quantity(3.0 * units.meter)
    Args = namedtuple("Args", "m r")
    return Args(m=m, r=r)


def test_basic_law(test_args):
    result = rotational_inertia_def.calculate_rotational_inertia(test_args.m, test_args.r)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.mass * units.length**2)
    result_value = convert_to(result, units.kilogram * units.meter**2).evalf(3)
    assert result_value == approx(9e-3, 1e-3)


def test_bad_mass(test_args):
    mb = Quantity(1.0 * units.second)
    with raises(errors.UnitsError):
        rotational_inertia_def.calculate_rotational_inertia(mb, test_args.r)
    with raises(TypeError):
        rotational_inertia_def.calculate_rotational_inertia(100, test_args.r)


def test_bad_radius(test_args):
    rb = Quantity(1.0 * units.second)
    with raises(errors.UnitsError):
        rotational_inertia_def.calculate_rotational_inertia(test_args.m, rb)
    with raises(TypeError):
        rotational_inertia_def.calculate_rotational_inertia(test_args.m, 100)
