from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.kinematic.rotational_inertia.geometries import (
    thin_rod_about_axis_through_center_perpendicular_to_length as thin_rod_formula)

# Description
## A thin rod of length of 30 cm and a mass of 0.5 kg rotates about an axis that passes through
## its center perpendicular to its length. Its rotaitonal inertia about that axis is 3.75e-3 kg*(m**2).

Args = namedtuple("Args", "m l")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    m = Quantity(0.5 * units.kilogram)
    l = Quantity(30.0 * units.centimeter)
    return Args(m=m, l=l)


def test_basic_law(test_args: Args) -> None:
    result = thin_rod_formula.calculate_rotational_inertia(test_args.m, test_args.l)
    assert_equal(result, 3.75e-3 * units.kilogram * units.meter**2)


def test_bad_mass(test_args: Args) -> None:
    mb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        thin_rod_formula.calculate_rotational_inertia(mb, test_args.l)
    with raises(TypeError):
        thin_rod_formula.calculate_rotational_inertia(100, test_args.l)


def test_bad_length(test_args: Args) -> None:
    lb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        thin_rod_formula.calculate_rotational_inertia(test_args.m, lb)
    with raises(TypeError):
        thin_rod_formula.calculate_rotational_inertia(test_args.m, 100)
