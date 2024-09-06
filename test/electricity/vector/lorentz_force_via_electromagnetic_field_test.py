from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal_vectors,
    assert_equal,
    units,
    prefixes,
    errors,
    Quantity,
    QuantityVector,
)
from symplyphysics.laws.electricity.vector import lorentz_force_via_electromagnetic_field as law

Args = namedtuple("Args", "f q e b v")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    f = QuantityVector([-7e-4 * units.newton, 2e-4 * units.newton, -5e-4 * units.newton])
    
    q = Quantity(-10 * prefixes.micro * units.coulomb)

    e_unit = units.volt / units.meter
    e = QuantityVector([100 * e_unit, 0, 100 * e_unit])

    v_unit = units.meter / units.second
    v = QuantityVector([-3 * v_unit, 2 * v_unit, 1 * v_unit])

    b = QuantityVector([10 * units.tesla, 10 * units.tesla, -10 * units.tesla])

    return Args(f=f, q=q, e=e, b=b, v=v)


def test_basic_law(test_args: Args) -> None:
    result = law.calculate_lorentz_force(test_args.q, test_args.e, test_args.b, test_args.v)
    assert_equal_vectors(result, test_args.f)


def test_electric_field_law(test_args: Args) -> None:
    result_vector = law.electric_field_law(
        lorentz_force_=test_args.f.to_base_vector(),
        magnetic_flux_density_=test_args.b.to_base_vector(),
        velocity_=test_args.v.to_base_vector(),
    )
    result = QuantityVector.from_base_vector(
        result_vector,
        subs={law.charge: test_args.q},
    )

    assert_equal_vectors(
        result,
        test_args.e,
        absolute_tolerance=1e-10,
    )


def test_charge_magnitude_law(test_args: Args) -> None:
    result = law.charge_law(
        lorentz_force_=test_args.f.to_base_vector(),
        electric_field_=test_args.e.to_base_vector(),
        magnetic_flux_density_=test_args.b.to_base_vector(),
        velocity_=test_args.v.to_base_vector(),
    )
    assert_equal(result, abs(test_args.q))


def test_bad_charge(test_args: Args) -> None:
    qb = Quantity(units.second)
    with raises(errors.UnitsError):
        law.calculate_lorentz_force(qb, test_args.e, test_args.b, test_args.v)
    with raises(TypeError):
        law.calculate_lorentz_force(100, test_args.e, test_args.b, test_args.v)


def test_bad_electric_field(test_args: Args) -> None:
    eb_vector = QuantityVector([units.second])
    with raises(errors.UnitsError):
        law.calculate_lorentz_force(test_args.q, eb_vector, test_args.b, test_args.v)

    eb_scalar = Quantity(units.volt / units.meter)
    with raises(AttributeError):
        law.calculate_lorentz_force(test_args.q, eb_scalar, test_args.b, test_args.v)

    with raises(TypeError):
        law.calculate_lorentz_force(test_args.q, 100, test_args.b, test_args.v)
    with raises(TypeError):
        law.calculate_lorentz_force(test_args.q, [100], test_args.b, test_args.v)


def test_bad_magnetic_field(test_args: Args) -> None:
    bb_vector = QuantityVector([units.second])
    with raises(errors.UnitsError):
        law.calculate_lorentz_force(test_args.q, test_args.e, bb_vector, test_args.v)

    bb_scalar = Quantity(units.tesla)
    with raises(AttributeError):
        law.calculate_lorentz_force(test_args.q, test_args.e, bb_scalar, test_args.v)

    with raises(TypeError):
        law.calculate_lorentz_force(test_args.q, test_args.e, 100, test_args.v)
    with raises(TypeError):
        law.calculate_lorentz_force(test_args.q, test_args.e, [100], test_args.v)


def test_bad_velocity(test_args: Args) -> None:
    vb_vector = QuantityVector([units.second])
    with raises(errors.UnitsError):
        law.calculate_lorentz_force(test_args.q, test_args.e, test_args.b, vb_vector)

    vb_scalar = Quantity(units.meter / units.second)
    with raises(AttributeError):
        law.calculate_lorentz_force(test_args.q, test_args.e, test_args.b, vb_scalar)

    with raises(TypeError):
        law.calculate_lorentz_force(test_args.q, test_args.e, test_args.b, 100)
    with raises(TypeError):
        law.calculate_lorentz_force(test_args.q, test_args.e, test_args.b, [100])
