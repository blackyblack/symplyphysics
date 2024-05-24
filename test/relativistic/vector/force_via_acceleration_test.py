from collections import namedtuple
from pytest import fixture, raises
from sympy.physics.units import speed_of_light
from symplyphysics import (
    assert_equal,
    units,
    Quantity,
    QuantityVector,
    errors,
)
from symplyphysics.laws.relativistic.vector import force_via_acceleration as law

Args = namedtuple("Args", "m a v f")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    m = Quantity(units.electron_rest_mass)
    a = QuantityVector([
        Quantity(1.0 * units.kilometer / units.second**2),
        Quantity(-2.0 * units.kilometer / units.second**2),
        Quantity(1.5 * units.kilometer / units.second**2),
    ])
    v = QuantityVector([
        Quantity(speed_of_light / 2),
        Quantity(speed_of_light / 2),
        Quantity(-1 * speed_of_light / 4),
    ])
    f = QuantityVector([
        Quantity(0),
        Quantity(-4.13e-27 * units.newton),
        Quantity(2.75e-27 * units.newton),
    ])
    return Args(m=m, a=a, v=v, f=f)


def test_force_law(test_args: Args) -> None:
    result = law.calculate_force(test_args.m, test_args.a, test_args.v)
    # TODO Finish
