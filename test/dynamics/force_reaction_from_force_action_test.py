from collections import namedtuple
from pytest import approx, fixture, raises
from sympy.vector import CoordSys3D

from symplyphysics import (
    solve, units, convert_to, SI, errors
)
from symplyphysics.laws.dynamics import force_reaction_from_force_action as newton_third_law

@fixture
def test_args():
    Faction = units.Quantity('Faction')
    SI.set_quantity_dimension(Faction, units.force)
    SI.set_quantity_scale_factor(Faction, 2 * units.newton)

    Args = namedtuple('Args', ['F'])
    return Args(F = Faction)

def test_basic_force():
    cartesian_coordinates = CoordSys3D('cartesian_coordinates')
    Faction = 2 * cartesian_coordinates.x + cartesian_coordinates.y

    result_force = solve(newton_third_law.law, newton_third_law.force_reaction, dict=True)[0][newton_third_law.force_reaction]
    result = result_force.subs(newton_third_law.force_action, Faction)

    # force action and force reaction should compensate each other
    assert (result + Faction) == 0
    # vector components should compensate each other
    assert result.coeff(cartesian_coordinates.x) == -1 * Faction.coeff(cartesian_coordinates.x)
    assert result.coeff(cartesian_coordinates.y) == -1 * Faction.coeff(cartesian_coordinates.y)
    assert result.coeff(cartesian_coordinates.z) == -1 * Faction.coeff(cartesian_coordinates.z)

def test_basic_force_quantity(test_args):
    result = newton_third_law.calculate_force_reaction(test_args.F)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.force)

    result_force = convert_to(result, units.newton).subs(units.newton, 1).evalf(2)
    assert result_force == approx(2.0, 0.01)

def test_bad_force():
    Fb = units.Quantity('Fb')
    SI.set_quantity_dimension(Fb, units.length)
    SI.set_quantity_scale_factor(Fb, 1 * units.meter)

    with raises(errors.UnitsError):
        newton_third_law.calculate_force_reaction(Fb)

    with raises(TypeError):
        newton_third_law.calculate_force_reaction(100)
