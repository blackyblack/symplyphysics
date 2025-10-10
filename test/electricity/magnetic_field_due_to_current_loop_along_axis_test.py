from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, units, Quantity, errors)
from symplyphysics.laws.electricity import magnetic_field_due_to_current_loop_along_axis as law

Args = namedtuple("Args", "i r z")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    i = Quantity(1 * units.ampere)
    r = Quantity(5e-2 * units.meter)
    z = Quantity(1 * units.meter)

    return Args(i=i, r=r, z=z)


def test_law(test_args: Args) -> None:
    result = law.calculate_magnetic_flux_density(test_args.i, test_args.r, test_args.z)
    assert_equal(result, 1.56e-9 * units.tesla, relative_tolerance=4e-3)


def test_bad_current(test_args: Args) -> None:
    ib = units.candela

    with raises(errors.UnitsError):
        law.calculate_magnetic_flux_density(ib, test_args.r, test_args.z)
    with raises(TypeError):
        law.calculate_magnetic_flux_density(100, test_args.r, test_args.z)


def test_bad_distance(test_args: Args) -> None:
    rb = units.candela

    with raises(errors.UnitsError):
        law.calculate_magnetic_flux_density(test_args.i, rb, test_args.z)
    with raises(TypeError):
        law.calculate_magnetic_flux_density(test_args.i, 100, test_args.z)

    with raises(errors.UnitsError):
        law.calculate_magnetic_flux_density(test_args.i, test_args.r, rb)
    with raises(TypeError):
        law.calculate_magnetic_flux_density(test_args.i, test_args.r, 100)
