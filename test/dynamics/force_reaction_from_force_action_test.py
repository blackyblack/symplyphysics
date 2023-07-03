from collections import namedtuple
from pytest import approx, fixture, raises
from sympy import solve
from sympy.vector import CoordSys3D
from symplyphysics import (
    errors,
    units,
    convert_to,
    Quantity,
    SI,
)
from symplyphysics.laws.dynamics import force_reaction_from_force_action as newton_third_law


@fixture(name="test_args")
def test_args_fixture():
    Fa = Quantity(2 * units.newton)
    Args = namedtuple("Args", ["F"])
    return Args(F=Fa)


def test_basic_force():
    cartesian_coordinates = CoordSys3D("cartesian_coordinates")
    # Make linter happy
    x = getattr(cartesian_coordinates, "x")
    y = getattr(cartesian_coordinates, "y")
    z = getattr(cartesian_coordinates, "z")
    Fa = 2 * x + y
    result_force = solve(newton_third_law.law, newton_third_law.force_reaction,
        dict=True)[0][newton_third_law.force_reaction]
    result = result_force.subs(newton_third_law.force_action, Fa)
    # force action and force reaction should compensate each other
    assert (result + Fa) == 0
    # vector components should compensate each other
    assert result.coeff(x) == -1 * Fa.coeff(x)
    assert result.coeff(y) == -1 * Fa.coeff(y)
    assert result.coeff(z) == -1 * Fa.coeff(z)


def test_basic_force_quantity(test_args):
    result = newton_third_law.calculate_force_reaction(test_args.F)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.force)
    result_force = convert_to(result, units.newton).evalf(2)
    assert result_force == approx(2.0, 0.01)


def test_bad_force():
    Fb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        newton_third_law.calculate_force_reaction(Fb)
    with raises(TypeError):
        newton_third_law.calculate_force_reaction(100)
