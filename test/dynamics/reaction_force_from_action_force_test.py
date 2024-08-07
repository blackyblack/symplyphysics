from collections import namedtuple
from pytest import fixture, raises
from sympy import solve
from sympy.vector import CoordSys3D
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.dynamics import reaction_force_from_action_force as newton_third_law

Args = namedtuple("Args", ["F"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    Fa = Quantity(2 * units.newton)
    return Args(F=Fa)


def test_basic_force() -> None:
    cartesian_coordinates = CoordSys3D("cartesian_coordinates")
    # Make linter happy
    x = getattr(cartesian_coordinates, "x")
    y = getattr(cartesian_coordinates, "y")
    z = getattr(cartesian_coordinates, "z")
    Fa = 2 * x + y
    result_force = solve(newton_third_law.law, newton_third_law.reaction_force,
        dict=True)[0][newton_third_law.reaction_force]
    result = result_force.subs(newton_third_law.action_force, Fa)
    # force action and force reaction should compensate each other
    assert (result + Fa) == 0
    # vector components should compensate each other
    assert result.coeff(x) == -1 * Fa.coeff(x)
    assert result.coeff(y) == -1 * Fa.coeff(y)
    assert result.coeff(z) == -1 * Fa.coeff(z)


def test_basic_force_quantity(test_args: Args) -> None:
    result = newton_third_law.calculate_force_reaction(test_args.F)
    assert_equal(result, 2 * units.newton)


def test_bad_force() -> None:
    Fb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        newton_third_law.calculate_force_reaction(Fb)
    with raises(TypeError):
        newton_third_law.calculate_force_reaction(100)
