from collections import namedtuple
from pytest import fixture, raises
from sympy.physics.units import speed_of_light, electron_rest_mass
from symplyphysics import (
    units,
    Quantity,
    QuantityVector,
    errors,
    assert_equal_vectors,
)
from symplyphysics.laws.relativistic.vector import relativistic_momentum as law

Args = namedtuple("Args", "m v p")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    m = Quantity(electron_rest_mass)
    v = QuantityVector([
        Quantity(5.0e-2 * speed_of_light),
        Quantity(1.0e-1 * speed_of_light),
        Quantity(-5.0e-1 * speed_of_light),
    ])
    p = QuantityVector([
        Quantity(5.82e-2 * electron_rest_mass * speed_of_light),
        Quantity(1.16e-1 * electron_rest_mass * speed_of_light),
        Quantity(-5.82e-1 * electron_rest_mass * speed_of_light),
    ])
    return Args(m=m, v=v, p=p)


def test_momentum_law(test_args: Args) -> None:
    result = law.calculate_momentum(test_args.m, test_args.v)
    assert_equal_vectors(result, test_args.p, relative_tolerance=4e-3)


def test_velocity_law(test_args: Args) -> None:
    result_vector = law.velocity_law(test_args.p.to_base_vector())
    result = QuantityVector.from_base_vector(
        result_vector,
        subs={law.rest_mass: test_args.m},
    )
    assert_equal_vectors(test_args.v, result, relative_tolerance=4e-3)


def test_bad_mass(test_args: Args) -> None:
    mb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        law.calculate_momentum(mb, test_args.v)
    with raises(TypeError):
        law.calculate_momentum(100, test_args.v)


def test_bad_velocity(test_args: Args) -> None:
    vb_vector = QuantityVector([Quantity(1 * units.coulomb)])
    with raises(errors.UnitsError):
        law.calculate_momentum(test_args.m, vb_vector)

    vb_scalar = Quantity(units.speed_of_light)
    with raises(AttributeError):
        law.calculate_momentum(test_args.m, vb_scalar)

    with raises(TypeError):
        law.calculate_momentum(test_args.m, 100)
    with raises(TypeError):
        law.calculate_momentum(test_args.m, [100])
