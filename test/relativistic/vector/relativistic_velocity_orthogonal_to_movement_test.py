from collections import namedtuple
from pytest import fixture, raises
from sympy.physics.units import speed_of_light
from symplyphysics import units, Quantity, errors
from symplyphysics.laws.relativistic.vector import relativistic_velocity_orthogonal_to_movement as law

from symplyphysics.core.experimental.vectors import VectorDot
from symplyphysics.core.experimental.coordinate_systems import CARTESIAN, QuantityCoordinateVector
from symplyphysics.core.experimental.approx import assert_equal_vectors

# Description
## The velocity vector of the body in the proper frame of reference is (0.1*c, 0.2*c, -0.333*c). The velocity vector of the
## proper frame of reference relative to the laboratory one is (-0.2*c, -0.4*c, 0). The velocity of the body in the lab frame
## of reference is (-0.111*c, -0.222*c, -0.331*c).

Args = namedtuple("Args", "up ul v")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    up = QuantityCoordinateVector([
        Quantity(0.100 * speed_of_light),
        Quantity(0.200 * speed_of_light),
        Quantity(-0.333 * speed_of_light)
    ], CARTESIAN)
    ul = QuantityCoordinateVector([
        Quantity(-0.111 * speed_of_light),
        Quantity(-0.222 * speed_of_light),
        Quantity(-0.331 * speed_of_light)
    ], CARTESIAN)
    v = QuantityCoordinateVector([
        Quantity(-0.200 * speed_of_light),
        Quantity(-0.400 * speed_of_light),
        Quantity(0),
    ], CARTESIAN)
    return Args(up=up, ul=ul, v=v)


def test_lab_law(test_args: Args) -> None:
    result = law.calculate_orthogonal_velocity_component_in_lab_frame(test_args.up, test_args.v)

    ul_tangential = (VectorDot(test_args.ul, test_args.v) / VectorDot(test_args.v, test_args.v) *
        test_args.v)
    ul_normal = test_args.ul - ul_tangential

    assert_equal_vectors(
        result,
        QuantityCoordinateVector.from_expr(ul_normal),
        relative_tolerance=3e-3,
        absolute_tolerance=1e-8,
    )


def test_bad_velocity(test_args: Args) -> None:
    vb_vector = QuantityCoordinateVector([Quantity(1 * units.coulomb), 0, 0], CARTESIAN)
    with raises(errors.UnitsError):
        law.calculate_orthogonal_velocity_component_in_lab_frame(vb_vector, test_args.v)
    with raises(errors.UnitsError):
        law.calculate_orthogonal_velocity_component_in_lab_frame(test_args.up, vb_vector)

    vb_scalar = Quantity(units.speed_of_light)
    with raises(ValueError):
        law.calculate_orthogonal_velocity_component_in_lab_frame(vb_scalar, test_args.v)
    with raises(ValueError):
        law.calculate_orthogonal_velocity_component_in_lab_frame(test_args.up, vb_scalar)

    with raises(TypeError):
        law.calculate_orthogonal_velocity_component_in_lab_frame(100, test_args.v)
    with raises(TypeError):
        law.calculate_orthogonal_velocity_component_in_lab_frame([100], test_args.v)
    with raises(TypeError):
        law.calculate_orthogonal_velocity_component_in_lab_frame(test_args.up, 100)
    with raises(TypeError):
        law.calculate_orthogonal_velocity_component_in_lab_frame(test_args.up, [100])
