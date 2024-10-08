from collections import namedtuple
from pytest import fixture, raises
from sympy import solve
from sympy.vector import CoordSys3D
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
    quantities,
)
from symplyphysics.laws.dynamics import buoyant_force_from_density_and_volume as archimedes_law

Args = namedtuple("Args", ["V", "pf"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    # water density
    pf = Quantity(1000 * units.kilogram / units.meter**3)
    V = Quantity(0.2 * units.meter**3)
    return Args(V=V, pf=pf)


def test_basic_force(test_args: Args) -> None:
    result = archimedes_law.calculate_force_buoyant(test_args.pf, test_args.V)
    assert_equal(result, 1961.3 * units.newton)


def test_force_vector() -> None:
    cartesian_coordinates = CoordSys3D("cartesian_coordinates")
    # Make linter happy
    x = getattr(cartesian_coordinates, "x")
    y = getattr(cartesian_coordinates, "y")
    z = getattr(cartesian_coordinates, "z")
    Fg = 2 * x + y
    result_force_expr = solve(archimedes_law.law, archimedes_law.buoyant_force,
        dict=True)[0][archimedes_law.buoyant_force]
    # set non vector variables to 1 since vector components cannot be properly calculated with quantities
    result = result_force_expr.subs({
        quantities.acceleration_due_to_gravity: Fg,
        archimedes_law.fluid_density: 1,
        archimedes_law.displaced_volume: 1
    })
    # make sure that result vector is not 0 (see Fg for non zero vector components)
    assert result.coeff(x) != 0
    assert result.coeff(y) != 0
    assert result.coeff(z) == 0
    # force action and force reaction should compensate each other
    assert (result + Fg) == 0
    # vector components should compensate each other (since density and volume = 1)
    assert result.coeff(x) == -1 * Fg.coeff(x)
    assert result.coeff(y) == -1 * Fg.coeff(y)
    assert result.coeff(z) == -1 * Fg.coeff(z)


def test_bad_density(test_args: Args) -> None:
    pb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        archimedes_law.calculate_force_buoyant(pb, test_args.V)
    with raises(TypeError):
        archimedes_law.calculate_force_buoyant(100, test_args.V)


def test_bad_volume(test_args: Args) -> None:
    Vb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        archimedes_law.calculate_force_buoyant(test_args.pf, Vb)
    with raises(TypeError):
        archimedes_law.calculate_force_buoyant(test_args.pf, 100)
