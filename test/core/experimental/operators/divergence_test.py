from sympy import sin, cos, ImmutableMatrix

from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.core.experimental.coordinate_systems import (CartesianCoordinateSystem,
    CylindricalCoordinateSystem, SphericalCoordinateSystem)
from symplyphysics.core.experimental.coordinate_systems.vector import CoordinateVector
from symplyphysics.core.experimental.coordinate_systems.scalar import CoordinateScalar
from symplyphysics.core.experimental.operators import (VectorGradient, VectorLaplacian,
    VectorDivergence)
from symplyphysics.core.experimental.solvers import vector_equals


def test_cartesian_system() -> None:
    cartesian = CartesianCoordinateSystem()

    x, y, z = cartesian.base_scalars

    vector = CoordinateVector([x, y, -z], cartesian)

    assert expr_equals(VectorDivergence(vector).scalar, 1)
