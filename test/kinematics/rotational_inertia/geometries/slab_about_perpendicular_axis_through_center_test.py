from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.kinematics.rotational_inertia.geometries import (
    slab_about_perpendicular_axis_through_center as slab_formula)

# Description
## A slab of length 40 cm, width of 30 cm and a mass of 5.0 kg is rotating an axis that passes
## through its center so that its length and width are perpendicular to it. Its rotational
## inertia amounts to about 0.1042 kg*(m**2).

Args = namedtuple("Args", "m a b")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    m = Quantity(5.0 * units.kilogram)
    a = Quantity(40.0 * units.centimeter)
    b = Quantity(30.0 * units.centimeter)
    return Args(m=m, a=a, b=b)


def test_basic_law(test_args: Args) -> None:
    result = slab_formula.calculate_rotational_inertia(test_args.m, test_args.a, test_args.b)
    assert_equal(result, 0.1042 * units.kilogram * units.meter**2)


def test_bad_mass(test_args: Args) -> None:
    mb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        slab_formula.calculate_rotational_inertia(mb, test_args.a, test_args.b)
    with raises(TypeError):
        slab_formula.calculate_rotational_inertia(100, test_args.a, test_args.b)


def test_bad_sizes(test_args: Args) -> None:
    lb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        slab_formula.calculate_rotational_inertia(test_args.m, lb, test_args.b)
    with raises(errors.UnitsError):
        slab_formula.calculate_rotational_inertia(test_args.m, test_args.a, lb)
    with raises(TypeError):
        slab_formula.calculate_rotational_inertia(test_args.m, 100, test_args.b)
    with raises(TypeError):
        slab_formula.calculate_rotational_inertia(test_args.m, test_args.a, 100)
