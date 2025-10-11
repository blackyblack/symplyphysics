from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import units, errors, Quantity, quantities
from symplyphysics.laws.electricity.vector import magnetic_field_due_to_constant_filamentary_current as law

from symplyphysics.core.experimental.coordinate_systems import CARTESIAN, QuantityCoordinateVector
from symplyphysics.core.experimental.approx import assert_equal_vectors

Args = namedtuple("Args", "mu i r l dl")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    mu = 1.1 * quantities.vacuum_permeability
    i = Quantity(1 * units.ampere)
    r = QuantityCoordinateVector([1 * units.meter, 0, 0], CARTESIAN)
    l = QuantityCoordinateVector([0, 0, 1 * units.meter], CARTESIAN)
    dl = QuantityCoordinateVector([0, 1e-3 * units.meter, 0], CARTESIAN)
    return Args(mu=mu, i=i, r=r, l=l, dl=dl)


def test_law(test_args: Args) -> None:
    result: QuantityCoordinateVector = law.calculate_magnetic_flux_density_change(
        test_args.mu, test_args.i, test_args.r, test_args.l, test_args.dl)
    expected = QuantityCoordinateVector([
        -3.89e-11 * units.tesla,
        0,
        -3.89e-11 * units.tesla,
    ], CARTESIAN)
    assert_equal_vectors(result, expected, relative_tolerance=2e-3)


def test_bad_permeability(test_args: Args) -> None:
    mub = units.meter
    with raises(errors.UnitsError):
        law.calculate_magnetic_flux_density_change(mub, test_args.i, test_args.r, test_args.l,
            test_args.dl)
    with raises(TypeError):
        law.calculate_magnetic_flux_density_change(100, test_args.i, test_args.r, test_args.l,
            test_args.dl)


def test_bad_current(test_args: Args) -> None:
    ib = units.meter
    with raises(errors.UnitsError):
        law.calculate_magnetic_flux_density_change(test_args.mu, ib, test_args.r, test_args.l,
            test_args.dl)
    with raises(TypeError):
        law.calculate_magnetic_flux_density_change(test_args.mu, 100, test_args.r, test_args.l,
            test_args.dl)


def test_bad_position_vector(test_args: Args) -> None:
    rb_vector = QuantityCoordinateVector([units.meter**2, 0, 0], CARTESIAN)
    with raises(errors.UnitsError):
        law.calculate_magnetic_flux_density_change(test_args.mu, test_args.i, rb_vector,
            test_args.l, test_args.dl)
    with raises(errors.UnitsError):
        law.calculate_magnetic_flux_density_change(test_args.mu, test_args.i, test_args.r,
            rb_vector, test_args.dl)
    with raises(errors.UnitsError):
        law.calculate_magnetic_flux_density_change(test_args.mu, test_args.i, test_args.r,
            test_args.l, rb_vector)

    rb_scalar = units.meter
    with raises(ValueError):
        law.calculate_magnetic_flux_density_change(test_args.mu, test_args.i, rb_scalar,
            test_args.l, test_args.dl)
    with raises(ValueError):
        law.calculate_magnetic_flux_density_change(test_args.mu, test_args.i, test_args.r,
            rb_scalar, test_args.dl)
    with raises(ValueError):
        law.calculate_magnetic_flux_density_change(test_args.mu, test_args.i, test_args.r,
            test_args.l, rb_scalar)

    with raises(TypeError):
        law.calculate_magnetic_flux_density_change(test_args.mu, test_args.i, 100, test_args.l,
            test_args.dl)
    with raises(TypeError):
        law.calculate_magnetic_flux_density_change(test_args.mu, test_args.i, [100], test_args.l,
            test_args.dl)
    with raises(TypeError):
        law.calculate_magnetic_flux_density_change(test_args.mu, test_args.i, test_args.r, 100,
            test_args.dl)
    with raises(TypeError):
        law.calculate_magnetic_flux_density_change(test_args.mu, test_args.i, test_args.r, [100],
            test_args.dl)
    with raises(TypeError):
        law.calculate_magnetic_flux_density_change(test_args.mu, test_args.i, test_args.r,
            test_args.l, 100)
    with raises(TypeError):
        law.calculate_magnetic_flux_density_change(test_args.mu, test_args.i, test_args.r,
            test_args.l, [100])
