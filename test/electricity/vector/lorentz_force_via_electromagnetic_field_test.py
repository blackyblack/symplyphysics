from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import units, prefixes, errors, Quantity
from symplyphysics.laws.electricity.vector import lorentz_force_via_electromagnetic_field as law

from symplyphysics.core.experimental.coordinate_systems import CARTESIAN, QuantityCoordinateVector
from symplyphysics.core.experimental.approx import assert_equal_vectors
from symplyphysics.core.experimental.solvers import solve_for_vector

Args = namedtuple("Args", "f q e b v")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    f = QuantityCoordinateVector([
        -7e-4 * units.newton,
        2e-4 * units.newton,
        -5e-4 * units.newton,
    ], CARTESIAN)

    q = Quantity(-10 * prefixes.micro * units.coulomb)

    e_unit = units.volt / units.meter
    e = QuantityCoordinateVector([
        100 * e_unit,
        0,
        100 * e_unit,
    ], CARTESIAN)

    v_unit = units.meter / units.second
    v = QuantityCoordinateVector([
        -3 * v_unit,
        2 * v_unit,
        1 * v_unit,
    ], CARTESIAN)

    b = QuantityCoordinateVector([
        10 * units.tesla,
        10 * units.tesla,
        -10 * units.tesla,
    ], CARTESIAN)

    return Args(f=f, q=q, e=e, b=b, v=v)


def test_basic_law(test_args: Args) -> None:
    result = law.calculate_lorentz_force(test_args.q, test_args.e, test_args.b, test_args.v)
    assert_equal_vectors(result, test_args.f)


def test_electric_field_law(test_args: Args) -> None:
    result = solve_for_vector(law.law, law.electric_field).subs({
        law.charge: test_args.q,
        law.lorentz_force: test_args.f,
        law.magnetic_flux_density: test_args.b,
        law.velocity: test_args.v,
    })

    result = QuantityCoordinateVector.from_expr(result)

    assert_equal_vectors(
        result,
        test_args.e,
        absolute_tolerance=1e-10,
    )


def test_bad_charge(test_args: Args) -> None:
    qb = Quantity(units.second)
    with raises(errors.UnitsError):
        law.calculate_lorentz_force(qb, test_args.e, test_args.b, test_args.v)
    with raises(TypeError):
        law.calculate_lorentz_force(100, test_args.e, test_args.b, test_args.v)


def test_bad_electric_field(test_args: Args) -> None:
    eb_vector = QuantityCoordinateVector([units.second, 0, 0], CARTESIAN)
    with raises(errors.UnitsError):
        law.calculate_lorentz_force(test_args.q, eb_vector, test_args.b, test_args.v)

    eb_scalar = Quantity(units.volt / units.meter)
    with raises(ValueError):
        law.calculate_lorentz_force(test_args.q, eb_scalar, test_args.b, test_args.v)

    with raises(TypeError):
        law.calculate_lorentz_force(test_args.q, 100, test_args.b, test_args.v)
    with raises(TypeError):
        law.calculate_lorentz_force(test_args.q, [100], test_args.b, test_args.v)


def test_bad_magnetic_field(test_args: Args) -> None:
    bb_vector = QuantityCoordinateVector([units.second, 0, 0], CARTESIAN)
    with raises(errors.UnitsError):
        law.calculate_lorentz_force(test_args.q, test_args.e, bb_vector, test_args.v)

    bb_scalar = Quantity(units.tesla)
    with raises(ValueError):
        law.calculate_lorentz_force(test_args.q, test_args.e, bb_scalar, test_args.v)

    with raises(TypeError):
        law.calculate_lorentz_force(test_args.q, test_args.e, 100, test_args.v)
    with raises(TypeError):
        law.calculate_lorentz_force(test_args.q, test_args.e, [100], test_args.v)


def test_bad_velocity(test_args: Args) -> None:
    vb_vector = QuantityCoordinateVector([units.second, 0, 0], CARTESIAN)
    with raises(errors.UnitsError):
        law.calculate_lorentz_force(test_args.q, test_args.e, test_args.b, vb_vector)

    vb_scalar = Quantity(units.meter / units.second)
    with raises(ValueError):
        law.calculate_lorentz_force(test_args.q, test_args.e, test_args.b, vb_scalar)

    with raises(TypeError):
        law.calculate_lorentz_force(test_args.q, test_args.e, test_args.b, 100)
    with raises(TypeError):
        law.calculate_lorentz_force(test_args.q, test_args.e, test_args.b, [100])
