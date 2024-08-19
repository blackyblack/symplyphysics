from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.kinematics.rotational_inertia import (
    rotational_inertia_about_axis_and_through_center_of_mass as parallel_axis_theorem,)

# Description
## A body weighing 3 kg is rotating about an axis. Its rotational inertia about a parallel
## axis that extends through the body's center of mass is 5.0 kg*(m**2). The two axis are 0.5 m apart.
## The rotational inertia about the original axis is 5.75 kg*(m**2).

Args = namedtuple("Args", "i_com m h")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    i_com = Quantity(5.0 * units.kilogram * units.meter**2)
    m = Quantity(3.0 * units.kilogram)
    h = Quantity(0.5 * units.meter)
    return Args(i_com=i_com, m=m, h=h)


def test_basic_law(test_args: Args) -> None:
    result = parallel_axis_theorem.calculate_rotational_inertia(test_args.i_com, test_args.m,
        test_args.h)
    assert_equal(result, 5.75 * units.kilogram * units.meter**2)


def test_bad_inertia(test_args: Args) -> None:
    ib = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        parallel_axis_theorem.calculate_rotational_inertia(ib, test_args.m, test_args.h)
    with raises(TypeError):
        parallel_axis_theorem.calculate_rotational_inertia(100, test_args.m, test_args.h)


def test_bad_mass(test_args: Args) -> None:
    mb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        parallel_axis_theorem.calculate_rotational_inertia(test_args.i_com, mb, test_args.h)
    with raises(TypeError):
        parallel_axis_theorem.calculate_rotational_inertia(test_args.i_com, 100, test_args.h)


def test_bad_distance(test_args: Args) -> None:
    hb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        parallel_axis_theorem.calculate_rotational_inertia(test_args.i_com, test_args.m, hb)
    with raises(TypeError):
        parallel_axis_theorem.calculate_rotational_inertia(test_args.i_com, test_args.m, 100)
