from collections import namedtuple
from pytest import approx, fixture, raises
from sympy.vector import CoordSys3D

from symplyphysics import (
    units, convert_to, SI, errors, solve
)
from symplyphysics.core.symbols.quantities import Quantity
from symplyphysics.laws.dynamics import buoyant_force_from_density_and_volume as archimedes_law

@fixture
def test_args():
    # water density
    pf = Quantity(1000 * units.kilogram / units.meter**3)
    V = Quantity(0.2 * units.meter**3)
    Args = namedtuple("Args", ["V", "pf"])
    return Args(V=V, pf=pf)

def test_basic_force(test_args):
    result = archimedes_law.calculate_force_buoyant(test_args.pf, test_args.V)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.force)
    result_force = convert_to(result, units.newton).subs(units.newton, 1).evalf(4)
    assert result_force == approx(1961.3, 0.01)

def test_force_vector():
    cartesian_coordinates = CoordSys3D("cartesian_coordinates")
    Fgravity = 2 * cartesian_coordinates.x + cartesian_coordinates.y
    result_force_expr = solve(archimedes_law.law, archimedes_law.force_buoyant, dict=True)[0][archimedes_law.force_buoyant]
    # set non vector variables to 1 since vector components cannot be properly calculated with quantities
    result = result_force_expr.subs({
        archimedes_law.units.acceleration_due_to_gravity: Fgravity,
        archimedes_law.fluid_density: 1,
        archimedes_law.displaced_volume: 1})
    # make sure that result vector is not 0 (see Fgravity for non zero vector components)
    assert result.coeff(cartesian_coordinates.x) != 0
    assert result.coeff(cartesian_coordinates.y) != 0
    assert result.coeff(cartesian_coordinates.z) == 0
    # force action and force reaction should compensate each other
    assert (result + Fgravity) == 0
    # vector components should compensate each other (since density and volume = 1)
    assert result.coeff(cartesian_coordinates.x) == -1 * Fgravity.coeff(cartesian_coordinates.x)
    assert result.coeff(cartesian_coordinates.y) == -1 * Fgravity.coeff(cartesian_coordinates.y)
    assert result.coeff(cartesian_coordinates.z) == -1 * Fgravity.coeff(cartesian_coordinates.z)

def test_bad_density(test_args):
    pb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        archimedes_law.calculate_force_buoyant(pb, test_args.V)
    with raises(TypeError):
        archimedes_law.calculate_force_buoyant(100, test_args.V)

def test_bad_volume(test_args):
    Vb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        archimedes_law.calculate_force_buoyant(test_args.pf, Vb)
    with raises(TypeError):
        archimedes_law.calculate_force_buoyant(test_args.pf, 100)
