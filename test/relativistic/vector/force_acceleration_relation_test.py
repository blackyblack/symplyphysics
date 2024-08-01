from collections import namedtuple
from pytest import fixture, raises
from sympy.physics.units import speed_of_light
from symplyphysics import (
    assert_equal,
    assert_equal_vectors,
    units,
    Quantity,
    QuantityVector,
    errors,
)
from symplyphysics.laws.relativistic.vector import force_acceleration_relation as law

Args = namedtuple("Args", "m a v f")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    m = Quantity(units.electron_rest_mass)
    a = QuantityVector([
        Quantity(1.0 * units.kilometer / units.second**2),
        Quantity(-2.0 * units.kilometer / units.second**2),
        Quantity(1.5 * units.kilometer / units.second**2),
    ])
    v = QuantityVector([
        Quantity(speed_of_light / 2),
        Quantity(speed_of_light / 2),
        Quantity(-1 * speed_of_light / 4),
    ])
    f = QuantityVector([
        Quantity(0),
        Quantity(-4.13e-27 * units.newton),
        Quantity(2.75e-27 * units.newton),
    ])
    return Args(m=m, a=a, v=v, f=f)


def test_acceleration_law(test_args: Args) -> None:
    result_vector = law.acceleration_law(
        test_args.f.to_base_vector(),
        test_args.v.to_base_vector(),
    )
    result = QuantityVector.from_base_vector(
        result_vector,
        subs={law.rest_mass: test_args.m},
    )
    assert_equal_vectors(result, test_args.a, tolerance=3e-3)


def test_force_law(test_args: Args) -> None:
    result = law.calculate_force(test_args.m, test_args.a, test_args.v)
    assert_equal_vectors(result, test_args.f, absolute_tolerance=1e-20)


def test_mass_law(test_args: Args) -> None:
    result = law.rest_mass_law(
        test_args.f.to_base_vector(),
        test_args.a.to_base_vector(),
        test_args.v.to_base_vector(),
    )
    assert_equal(result, test_args.m, tolerance=1e-3)


def test_bad_mass(test_args: Args) -> None:
    mb = Quantity(units.coulomb)
    with raises(errors.UnitsError):
        law.calculate_force(mb, test_args.a, test_args.v)
    with raises(TypeError):
        law.calculate_force(100, test_args.a, test_args.v)


def test_bad_acceleration(test_args: Args) -> None:
    ab_vector = QuantityVector([Quantity(units.coulomb)])
    with raises(errors.UnitsError):
        law.calculate_force(test_args.m, ab_vector, test_args.v)

    ab_scalar = Quantity(units.planck_acceleration)
    with raises(AttributeError):
        law.calculate_force(test_args.m, ab_scalar, test_args.v)

    with raises(TypeError):
        law.calculate_force(test_args.m, 100, test_args.v)
    with raises(TypeError):
        law.calculate_force(test_args.m, [100], test_args.v)


def test_bad_velocity(test_args: Args) -> None:
    vb_vector = QuantityVector([Quantity(units.coulomb)])
    with raises(errors.UnitsError):
        law.calculate_force(test_args.m, test_args.a, vb_vector)

    vb_scalar = Quantity(units.speed_of_light)
    with raises(AttributeError):
        law.calculate_force(test_args.m, test_args.a, vb_scalar)

    with raises(TypeError):
        law.calculate_force(test_args.m, test_args.a, 100)
    with raises(TypeError):
        law.calculate_force(test_args.m, test_args.a, [100])
