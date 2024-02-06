from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_approx,
    errors,
    units,
    convert_to,
    Quantity,
    SI,
)
from symplyphysics.laws.kinematic.rotational_inertia.geometries import (
    thin_rod_about_axis_through_center_perpendicular_to_length as thin_rod_formula)

# Description
## A thin rod of length of 30 cm and a mass of 0.5 kg rotates about an axis that passes through
## its center perpendicular to its length. Its rotaitonal inertia about that axis is 3.75e-3 kg*(m**2).


@fixture(name="test_args")
def test_args_fixture():
    m = Quantity(0.5 * units.kilogram)
    l = Quantity(30.0 * units.centimeter)
    Args = namedtuple("Args", "m l")
    return Args(m=m, l=l)


def test_basic_law(test_args):
    result = thin_rod_formula.calculate_rotational_inertia(test_args.m, test_args.l)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.mass * units.length**2)
    result_value = convert_to(result, units.kilogram * units.meter**2).evalf(3)
    assert_approx(result_value, 3.75e-3)


def test_bad_mass(test_args):
    mb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        thin_rod_formula.calculate_rotational_inertia(mb, test_args.l)
    with raises(TypeError):
        thin_rod_formula.calculate_rotational_inertia(100, test_args.l)


def test_bad_length(test_args):
    lb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        thin_rod_formula.calculate_rotational_inertia(test_args.m, lb)
    with raises(TypeError):
        thin_rod_formula.calculate_rotational_inertia(test_args.m, 100)
