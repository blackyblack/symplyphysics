from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.hydro import shear_stress_is_shear_force_over_area as shear_stress

# The shear stress of a force of 1 N applied to an area of 1 m^2 is 1 Pa.


@fixture(name="test_args")
def test_args_fixture():
    F = Quantity(1 * units.newton)
    A = Quantity(1 * units.meter**2)
    Args = namedtuple("Args", "F A")
    return Args(F=F, A=A)


def test_basic_law(test_args):
    result = shear_stress.calculate_shear_stress(test_args.F, test_args.A)
    assert_equal(result, 1 * units.pascal)


def test_bad_force(test_args):
    Fb = Quantity(1 * units.second)
    with raises(errors.UnitsError):
        shear_stress.calculate_shear_stress(Fb, test_args.A)
    with raises(TypeError):
        shear_stress.calculate_shear_stress(100, test_args.A)


def test_bad_area(test_args):
    Ab = Quantity(1 * units.second)
    with raises(errors.UnitsError):
        shear_stress.calculate_shear_stress(test_args.F, Ab)
    with raises(TypeError):
        shear_stress.calculate_shear_stress(test_args.F, 100)
