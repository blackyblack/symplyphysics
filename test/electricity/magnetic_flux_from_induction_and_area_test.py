from collections import namedtuple
from pytest import fixture, raises
from sympy import pi
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.electricity import magnetic_flux_from_induction_and_area as flux_law

# Description
## With a magnetic induction of 0.3 tesla and an area of 0.1 meters, the magnetic flux is 2.12e-2 weber.
## The angle between the normal of the pad and the magnetic induction is 45 degree (pi / 4 radian).
## https://physics.icalculator.com/magnetic-flux-calculator.html


@fixture(name="test_args")
def test_args_fixture():
    induction = Quantity(0.3 * units.tesla)
    area = Quantity(0.1 * units.meter**2)
    angle = pi / 4

    Args = namedtuple("Args", ["induction", "area", "angle"])
    return Args(induction=induction, area=area, angle=angle)


def test_basic_flux(test_args):
    result = flux_law.calculate_flux(test_args.induction, test_args.area, test_args.angle)
    assert_equal(result, 2.12e-2 * units.weber)


def test_bad_induction(test_args):
    induction = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        flux_law.calculate_flux(induction, test_args.area, test_args.angle)
    with raises(TypeError):
        flux_law.calculate_flux(100, test_args.area, test_args.angle)


def test_bad_area(test_args):
    area = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        flux_law.calculate_flux(test_args.induction, area, test_args.angle)
    with raises(TypeError):
        flux_law.calculate_flux(test_args.induction, 100, test_args.angle)


def test_bad_angle(test_args):
    angle = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        flux_law.calculate_flux(test_args.induction, test_args.area, angle)
    with raises(AttributeError):
        flux_law.calculate_flux(test_args.induction, test_args.area, True)
