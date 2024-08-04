from collections import namedtuple
from pytest import fixture, raises
from sympy import Rational
from symplyphysics import (assert_equal, units, Quantity, errors, prefixes)
from symplyphysics.laws.electricity.circuits.diodes import diode_constant_of_cylindrical_diode as constant_law

# Description
## The area of the anode is 30 [centimeter^2]. The radius of the anode is 10 centimeter, the radius of the cathode
## is 0.5 centimeter. Then the diode constant is 0.776 [microampere / volt^(3 / 2)].

Args = namedtuple("Args", ["anode_area", "anode_radius", "cathode_radius"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    anode_area = Quantity(30 * units.centimeter**2)
    anode_radius = Quantity(10 * units.centimeter)
    cathode_radius = Quantity(0.5 * units.centimeter)

    return Args(anode_area=anode_area, anode_radius=anode_radius, cathode_radius=cathode_radius)


def test_basic_diode_constant(test_args: Args) -> None:
    result = constant_law.calculate_diode_constant(test_args.anode_area, test_args.anode_radius,
        test_args.cathode_radius)
    assert_equal(result, 0.776 * prefixes.micro * units.ampere / units.volt**Rational(3, 2))


def test_bad_anode_area(test_args: Args) -> None:
    anode_area = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        constant_law.calculate_diode_constant(anode_area, test_args.anode_radius,
            test_args.cathode_radius)
    with raises(TypeError):
        constant_law.calculate_diode_constant(100, test_args.anode_radius, test_args.cathode_radius)


def test_bad_anode_radius(test_args: Args) -> None:
    anode_radius = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        constant_law.calculate_diode_constant(test_args.anode_area, anode_radius,
            test_args.cathode_radius)
    with raises(TypeError):
        constant_law.calculate_diode_constant(test_args.anode_area, 100, test_args.cathode_radius)
    with raises(ValueError):
        constant_law.calculate_diode_constant(test_args.anode_area, test_args.cathode_radius,
            test_args.anode_radius)


def test_bad_cathode_radius(test_args: Args) -> None:
    cathode_radius = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        constant_law.calculate_diode_constant(test_args.anode_area, test_args.anode_radius,
            cathode_radius)
    with raises(TypeError):
        constant_law.calculate_diode_constant(test_args.anode_area, test_args.anode_radius, 100)
