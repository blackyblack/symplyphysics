from collections import namedtuple
from pytest import approx, fixture, raises
from sympy import pi
from symplyphysics import (units, SI, convert_to, Quantity, errors)
from symplyphysics.laws.electricity import ampere_law

# Description
## With a current of 3.5 amperes, a conductor length of 2 meters and a magnetic induction of 0.6 tesla,
## the force will be equal to 3 newton. The angle between the magnetic induction and the current direction
## is 45 degrees (pi / 4 radians).
## https://physics.icalculator.com/calculating-magnetic-field-using-the-amperes-law.html


@fixture(name="test_args")
def test_args_fixture():
    current = Quantity(3.5 * units.ampere)
    length = Quantity(2 * units.meter)
    induction = Quantity(0.6 * units.tesla)
    angle = pi/4

    Args = namedtuple("Args", ["current", "length", "induction", "angle"])
    return Args(current=current, length=length, angle=angle, induction=induction)


def test_basic_force(test_args):
    result = ampere_law.calculate_force(test_args.current, test_args.length, test_args.angle, test_args.induction)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.force)
    result = convert_to(result, units.newton).evalf(5)
    assert result == approx(3, rel=0.1)


def test_bad_current(test_args):
    current = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        ampere_law.calculate_force(current, test_args.length, test_args.angle, test_args.induction)
    with raises(TypeError):
        ampere_law.calculate_force(100, test_args.length, test_args.angle, test_args.induction)


def test_bad_length(test_args):
    length = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        ampere_law.calculate_force(test_args.current, length, test_args.angle, test_args.induction)
    with raises(TypeError):
        ampere_law.calculate_force(test_args.current, 100, test_args.angle, test_args.induction)


def test_bad_angle(test_args):
    angle = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        ampere_law.calculate_force(test_args.current, test_args.length, angle, test_args.induction)
    with raises(AttributeError):
        ampere_law.calculate_force(test_args.current, test_args.length, True, test_args.induction)


def test_bad_induction(test_args):
    induction = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        ampere_law.calculate_force(test_args.current, test_args.length, test_args.angle, induction)
    with raises(TypeError):
        ampere_law.calculate_force(test_args.current, test_args.length, test_args.angle, 100)
