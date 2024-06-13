from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    errors,
    units,
    Quantity,
    QuantityVector,
)
from symplyphysics.core.approx import assert_equal_vectors
from symplyphysics.laws.kinematic.vector import velocity_of_transfer_between_reference_frames as law

Args = namedtuple("Args", "vtr v0 w r")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    velocity_unit = units.meter / units.second
    vtr = QuantityVector([
        3 * velocity_unit,
        -3 * velocity_unit,
        6 * velocity_unit,
    ])
    v0 = QuantityVector([
        1 * velocity_unit,
        -1 * velocity_unit,
        4 * velocity_unit,
    ])
    w = QuantityVector([
        Quantity(1 * units.radian / units.second),
        Quantity(0),
        Quantity(-1 * units.radian / units.second),
    ])
    r = QuantityVector([
        3 * units.meter,
        2 * units.meter,
        -1 * units.meter,
    ])
    return Args(vtr=vtr, v0=v0, w=w, r=r)


def test_transfer_law(test_args: Args) -> None:
    result = law.calculate_transfer_velocity(test_args.v0, test_args.w, test_args.r)
    assert_equal_vectors(result, test_args.vtr)


def test_origin_law(test_args: Args) -> None:
    result_vector = law.origin_velocity_law(
        test_args.vtr.to_base_vector(),
        test_args.w.to_base_vector(),
        test_args.r.to_base_vector(),
    )
    result = QuantityVector.from_base_vector(result_vector)
    assert_equal_vectors(result, test_args.v0)


def test_bad_velocity(test_args: Args) -> None:
    vb_vector = QuantityVector([Quantity(1 * units.coulomb)])
    with raises(errors.UnitsError):
        law.calculate_transfer_velocity(vb_vector, test_args.w, test_args.r)

    vb_scalar = Quantity(1 * units.meter / units.second)
    with raises(AttributeError):
        law.calculate_transfer_velocity(vb_scalar, test_args.w, test_args.r)

    with raises(TypeError):
        law.calculate_transfer_velocity(100, test_args.w, test_args.r)
    with raises(TypeError):
        law.calculate_transfer_velocity([100], test_args.w, test_args.r)


def test_bad_angular_velocity(test_args: Args) -> None:
    wb_vector = QuantityVector([Quantity(1 * units.coulomb)])
    with raises(errors.UnitsError):
        law.calculate_transfer_velocity(test_args.v0, wb_vector, test_args.r)
    
    wb_scalar = Quantity(1 * units.radian / units.second)
    with raises(AttributeError):
        law.calculate_transfer_velocity(test_args.v0, wb_scalar, test_args.r)
    
    with raises(TypeError):
        law.calculate_transfer_velocity(test_args.v0, 100, test_args.r)
    with raises(TypeError):
        law.calculate_transfer_velocity(test_args.v0, [100], test_args.r)


def test_bad_position(test_args: Args) -> None:
    rb_vector = QuantityVector([Quantity(1 * units.coulomb)])
    with raises(errors.UnitsError):
        law.calculate_transfer_velocity(test_args.v0, test_args.w, rb_vector)
    
    rb_scalar = Quantity(1 * units.meter)
    with raises(AttributeError):
        law.calculate_transfer_velocity(test_args.v0, test_args.w, rb_scalar)
    
    with raises(TypeError):
        law.calculate_transfer_velocity(test_args.v0, test_args.w, 100)
    with raises(TypeError):
        law.calculate_transfer_velocity(test_args.v0, test_args.w, [100])
