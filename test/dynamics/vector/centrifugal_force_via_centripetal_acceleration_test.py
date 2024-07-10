from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    assert_equal_vectors,
    units,
    errors,
    Quantity,
    QuantityVector,
)
from symplyphysics.laws.dynamics.vector import centrifugal_force_via_centripetal_acceleration as law

Args = namedtuple("Args", "m a f")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    m = Quantity(2.0 * units.kilogram)
    
    a_unit = units.meter / units.second**2
    a = QuantityVector([0.5 * a_unit, -3.0 * a_unit, 4.5 * a_unit])

    f = QuantityVector([-1.0 * units.newton, 6.0 * units.newton, -9.0 * units.newton])

    return Args(m=m, a=a, f=f)


def test_force_law(test_args: Args) -> None:
    result = law.calculate_centrifugal_force(test_args.m, test_args.a)
    assert_equal_vectors(result, test_args.f)


def test_acceleration_law(test_args: Args) -> None:
    result_vector = law.centripetal_acceleration_law(test_args.f.to_base_vector())
    result = QuantityVector.from_base_vector(result_vector, subs={law.mass: test_args.m})
    assert_equal_vectors(result, test_args.a)


def test_mass_law(test_args: Args) -> None:
    result = law.mass_law(
        test_args.f.to_base_vector(),
        test_args.a.to_base_vector(),
    )
    assert_equal(result, test_args.m)


def test_bad_mass(test_args: Args) -> None:
    mb = Quantity(units.coulomb)
    with raises(errors.UnitsError):
        law.calculate_centrifugal_force(mb, test_args.a)
    with raises(TypeError):
        law.calculate_centrifugal_force(100, test_args.a)


def test_bad_acceleration(test_args: Args) -> None:
    ab_vector = QuantityVector([units.coulomb])
    with raises(errors.UnitsError):
        law.calculate_centrifugal_force(test_args.m, ab_vector)
    
    ab_scalar = units.meter / units.second**2
    with raises(AttributeError):
        law.calculate_centrifugal_force(test_args.m, ab_scalar)

    with raises(TypeError):
        law.calculate_centrifugal_force(test_args.m, 100)
    with raises(TypeError):
        law.calculate_centrifugal_force(test_args.m, [100])
