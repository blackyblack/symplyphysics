from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    errors,
    units,
    Quantity,
    QuantityVector,
    assert_equal_vectors,
)
from symplyphysics.laws.gravity.vector import falling_body_displacement as law

Args = namedtuple("Args", "t v0 w g")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    t = Quantity(5 * units.second)

    v_unit = units.meter / units.second
    v0 = QuantityVector([v_unit, v_unit, v_unit])

    w = QuantityVector([0, -7.27e-5 * units.radian / units.second, 0])

    g = QuantityVector([0, 0, -9.81 * units.meter / units.second**2])

    return Args(t=t, v0=v0, w=w, g=g)


def test_law(test_args: Args) -> None:
    result = law.calculate_displacement(test_args.t, test_args.v0, test_args.w, test_args.g)
    assert_equal_vectors(
        result,
        QuantityVector([4.97 * units.meter, 5.00 * units.meter, -118.0 * units.meter]),
        tolerance=4e-3,
    )


def test_bad_time(test_args: Args) -> None:
    tb = Quantity(units.coulomb)
    with raises(errors.UnitsError):
        law.calculate_displacement(tb, test_args.v0, test_args.w, test_args.g)
    with raises(TypeError):
        law.calculate_displacement(100, test_args.v0, test_args.w, test_args.g)


def test_bad_linear_velocity(test_args: Args) -> None:
    vb_vector = QuantityVector([units.coulomb])
    with raises(errors.UnitsError):
        law.calculate_displacement(test_args.t, vb_vector, test_args.w, test_args.g)

    vb_scalar = Quantity(units.meter / units.second)
    with raises(AttributeError):
        law.calculate_displacement(test_args.t, vb_scalar, test_args.w, test_args.g)

    with raises(TypeError):
        law.calculate_displacement(test_args.t, 100, test_args.w, test_args.g)
    with raises(TypeError):
        law.calculate_displacement(test_args.t, [100], test_args.w, test_args.g)


def test_bad_angular_velocity(test_args: Args) -> None:
    wb_vector = QuantityVector([units.coulomb])
    with raises(errors.UnitsError):
        law.calculate_displacement(test_args.t, test_args.v0, wb_vector, test_args.g)

    wb_scalar = Quantity(units.radian / units.second)
    with raises(AttributeError):
        law.calculate_displacement(test_args.t, test_args.v0, wb_scalar, test_args.g)

    with raises(TypeError):
        law.calculate_displacement(test_args.t, test_args.v0, 100, test_args.g)
    with raises(TypeError):
        law.calculate_displacement(test_args.t, test_args.v0, [100], test_args.g)


def test_bad_acceleration(test_args: Args) -> None:
    gb_vector = QuantityVector([units.coulomb])
    with raises(errors.UnitsError):
        law.calculate_displacement(test_args.t, test_args.v0, test_args.w, gb_vector)

    gb_scalar = Quantity(units.meter / units.second**2)
    with raises(AttributeError):
        law.calculate_displacement(test_args.t, test_args.v0, test_args.w, gb_scalar)

    with raises(TypeError):
        law.calculate_displacement(test_args.t, test_args.v0, test_args.w, 100)
    with raises(TypeError):
        law.calculate_displacement(test_args.t, test_args.v0, test_args.w, [100])
