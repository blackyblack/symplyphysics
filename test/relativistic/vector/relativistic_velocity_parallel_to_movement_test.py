from collections import namedtuple
from pytest import fixture, raises
from sympy.physics.units import speed_of_light
from symplyphysics import units, Quantity, errors
from symplyphysics.laws.relativistic.vector import relativistic_velocity_parallel_to_movement as law

from symplyphysics.core.experimental.vectors import VectorDot
from symplyphysics.core.experimental.coordinate_systems import CARTESIAN, QuantityCoordinateVector
from symplyphysics.core.experimental.approx import assert_equal_vectors

# Description
## A body is moving parallel to the velocity vector of the proper reference frame relative to the lab frame.
## The velocity vector of the body in its proper frame of reference is (0.1*c, 0.2*c, 0). The velocity vector of the
## proper frame of reference relative to the laboratory one is (-0.2*c, -0.4*c, 0). Then the velocity of the body relative
## to the laboratory frame of reference is (-0.111*c, -0.222*c, 0).

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
    result = law.calculate_parallel_velocity_component_in_lab_frame(test_args.up, test_args.v)

    ul_tangential = (VectorDot(test_args.ul, test_args.v) / VectorDot(test_args.v, test_args.v) *
        test_args.v)
    ul_tangential = QuantityCoordinateVector.from_expr(ul_tangential)

    assert_equal_vectors(result, ul_tangential, relative_tolerance=2e-3)


def test_proper_law(test_args: Args) -> None:
    v = QuantityCoordinateVector.from_expr(-test_args.v)
    result = law.calculate_parallel_velocity_component_in_lab_frame(test_args.ul, v)

    up_tangential = (VectorDot(test_args.up, test_args.v) / VectorDot(test_args.v, test_args.v) *
        test_args.v)
    up_tangential = QuantityCoordinateVector.from_expr(up_tangential)

    assert_equal_vectors(result, up_tangential, relative_tolerance=2e-3)


def test_bad_object_velocity(test_args: Args) -> None:
    qb = Quantity(1 * units.coulomb)
    vb_vector = QuantityCoordinateVector([qb, qb, qb], CARTESIAN)
    with raises(ValueError):
        law.calculate_parallel_velocity_component_in_lab_frame(vb_vector, test_args.v)

    vb_scalar = Quantity(units.speed_of_light)
    with raises(ValueError):
        law.calculate_parallel_velocity_component_in_lab_frame(vb_scalar, test_args.v)

    with raises(ValueError):
        law.calculate_parallel_velocity_component_in_lab_frame(100, test_args.v)
    with raises(ValueError):
        law.calculate_parallel_velocity_component_in_lab_frame([100], test_args.v)


def test_bad_frame_velocity(test_args: Args) -> None:
    qb = Quantity(1 * units.coulomb)
    vb_vector = QuantityCoordinateVector([qb, qb, qb], CARTESIAN)
    with raises(errors.UnitsError):
        law.calculate_parallel_velocity_component_in_lab_frame(test_args.up, vb_vector)

    vb_scalar = Quantity(units.speed_of_light)
    with raises(ValueError):
        law.calculate_parallel_velocity_component_in_lab_frame(test_args.up, vb_scalar)

    with raises(TypeError):
        law.calculate_parallel_velocity_component_in_lab_frame(test_args.up, 100)
    with raises(TypeError):
        law.calculate_parallel_velocity_component_in_lab_frame(test_args.up, [100])
