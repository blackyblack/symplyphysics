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
from symplyphysics.laws.relativistic.vector import relativistic_mass_moment as law

Args = namedtuple("Args", "m x v t")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    m = Quantity(1.0 * electron_rest_mass)
    x = QuantityVector([
        Quantity(100.0 * units.meter),
        Quantity(-10.0 * units.meter),
        Quantity(40.0 * units.meter),
    ])
    v = QuantityVector([
        Quantity(0.0001 * speed_of_light),
        Quantity(0.0002 * speed_of_light),
        Quantity(-0.0003 * speed_of_light),
    ])
    t = Quantity(1.0 * units.second)
    return Args(m=m, x=x, v=v, t=t)


def test_law(test_args: Args) -> None:
    result = law.calculate_mass_moment(test_args.m, test_args.x, test_args.v, test_args.t)

    measurement_unit = units.electron_rest_mass * units.kilometer
    correct = QuantityVector([
        Quantity(-29.9 * measurement_unit),
        Quantity(-60.0 * measurement_unit),
        Quantity(90.0 * measurement_unit),
    ])

    assert_equal_vectors(result, correct)


def test_bad_mass(test_args: Args) -> None:
    mb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        law.calculate_mass_moment(mb, test_args.x, test_args.v, test_args.t)
    with raises(TypeError):
        law.calculate_mass_moment(100, test_args.x, test_args.v, test_args.t)


def test_bad_position(test_args: Args) -> None:
    xb_vector = QuantityVector([Quantity(1 * units.coulomb)])
    with raises(errors.UnitsError):
        law.calculate_mass_moment(test_args.m, xb_vector, test_args.v, test_args.t)

    xb_scalar = Quantity(1 * units.meter)
    with raises(AttributeError):
        law.calculate_mass_moment(test_args.m, xb_scalar, test_args.v, test_args.t)

    with raises(TypeError):
        law.calculate_mass_moment(test_args.m, 100, test_args.v, test_args.t)
    with raises(TypeError):
        law.calculate_mass_moment(test_args.m, [100], test_args.v, test_args.t)


def test_bad_velocity(test_args: Args) -> None:
    vb_vector = QuantityVector([Quantity(1 * units.coulomb)])
    with raises(errors.UnitsError):
        law.calculate_mass_moment(test_args.m, test_args.x, vb_vector, test_args.t)

    vb_scalar = Quantity(1 * units.meter / units.second)
    with raises(AttributeError):
        law.calculate_mass_moment(test_args.m, test_args.x, vb_scalar, test_args.t)

    with raises(TypeError):
        law.calculate_mass_moment(test_args.m, test_args.x, 100, test_args.t)
    with raises(TypeError):
        law.calculate_mass_moment(test_args.m, test_args.x, [100], test_args.t)


def test_bad_time(test_args: Args) -> None:
    tb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        law.calculate_mass_moment(test_args.m, test_args.x, test_args.v, tb)
    with raises(TypeError):
        law.calculate_mass_moment(test_args.m, test_args.x, test_args.v, 100)
