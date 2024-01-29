from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    errors,
    units,
    convert_to,
    Quantity,
    SI,
)
from symplyphysics.laws.kinematic import (
    rotational_inertia_about_axis_and_through_com as parallel_axis_theorem,
)

# Description
## A body weighing 3 kg is rotating about an axis. Its rotational inertia about a parallel
## axis that extends through the body's com is 5.0 kg*(m**2). The two axis are 0.5 m apart.
## The rotational inertia about the original axis is 5.75 kg*(m**2).


@fixture(name="test_args")
def test_args_fixture():
    i_com = Quantity(5.0 * units.kilogram * units.meter**2)
    m = Quantity(3.0 * units.kilogram)
    h = Quantity(0.5 * units.meter)
    Args = namedtuple("Args", "i_com m h")
    return Args(i_com=i_com, m=m, h=h)


def test_basic_law(test_args):
    result = parallel_axis_theorem.calculate_rotational_inertia(
        test_args.i_com, test_args.m, test_args.h
    )
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.mass * units.length**2)
    result_value = convert_to(result, units.kilogram * units.meter**2).evalf(3)
    assert result_value == approx(5.75, 1e-3)


def test_bad_inertia(test_args):
    ib = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        parallel_axis_theorem.calculate_rotational_inertia(ib, test_args.m, test_args.h)
    with raises(TypeError):
        parallel_axis_theorem.calculate_rotational_inertia(100, test_args.m, test_args.h)


def test_bad_mass(test_args):
    mb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        parallel_axis_theorem.calculate_rotational_inertia(test_args.i_com, mb, test_args.h)
    with raises(TypeError):
        parallel_axis_theorem.calculate_rotational_inertia(test_args.i_com, 100, test_args.h)

def test_bad_distance(test_args):
    hb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        parallel_axis_theorem.calculate_rotational_inertia(test_args.i_com, test_args.m, hb)
    with raises(TypeError):
        parallel_axis_theorem.calculate_rotational_inertia(test_args.i_com, test_args.m, 100)
