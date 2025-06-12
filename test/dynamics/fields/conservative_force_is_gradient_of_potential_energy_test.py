from collections import namedtuple
from pytest import fixture
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.dynamics.fields import (
    conservative_force_is_gradient_of_potential_energy as gradient_law,)

from symplyphysics.core.experimental.coordinate_systems import (CARTESIAN, CoordinateVector,
    CoordinateScalar)
from symplyphysics.core.experimental.solvers import vector_equals

# Description
## The force associated with a potential of x^2/2 has the vector form of F = -x*e_x, where
## e_x is the unit vector in the direction of the x-axis

Args = namedtuple("Args", "potential")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    x, _, _ = CARTESIAN.base_scalars

    potential = CoordinateScalar(x**2 / 2, CARTESIAN)

    return Args(potential=potential)


def test_basic_law(test_args: Args) -> None:
    result = gradient_law.law.rhs.subs(
        gradient_law.potential_energy(gradient_law.position_vector),
        test_args.potential,
    ).doit()
    result = CoordinateVector.from_expr(result)

    x, _, _ = CARTESIAN.base_scalars
    expected = CoordinateVector([-x, 0, 0], CARTESIAN)

    assert vector_equals(result, expected)
