from collections import namedtuple
from pytest import fixture, raises
from sympy import pi
from symplyphysics import (assert_equal, units, Quantity, errors)
from symplyphysics.laws.electricity import magnetic_field_due_to_finite_coil_along_axis as law

Args = namedtuple("Args", "i n l p1 p2")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    i = Quantity(10 * units.ampere)
    n = 1000
    l = Quantity(1 * units.meter)
    p1 = pi / 6
    p2 = pi / 3
    return Args(i=i, n=n, l=l, p1=p1, p2=p2)


def test_law(test_args: Args) -> None:
    result = law.calculate_magnetic_flux_density(test_args.i, test_args.n, test_args.l,
        test_args.p1, test_args.p2)
    assert_equal(result, 8.6e-3 * units.tesla, relative_tolerance=2e-3)


def test_bad_current(test_args: Args) -> None:
    ib = units.candela

    with raises(errors.UnitsError):
        law.calculate_magnetic_flux_density(ib, test_args.n, test_args.l, test_args.p1,
            test_args.p2)
    with raises(TypeError):
        law.calculate_magnetic_flux_density(100, test_args.n, test_args.l, test_args.p1,
            test_args.p2)


def test_bad_length(test_args: Args) -> None:
    lb = units.candela

    with raises(errors.UnitsError):
        law.calculate_magnetic_flux_density(test_args.i, test_args.n, lb, test_args.p1,
            test_args.p2)
    with raises(TypeError):
        law.calculate_magnetic_flux_density(test_args.i, test_args.n, 100, test_args.p1,
            test_args.p2)


def test_bad_dimensionless(test_args: Args) -> None:
    nb = units.candela

    with raises(errors.UnitsError):
        law.calculate_magnetic_flux_density(test_args.i, nb, test_args.l, test_args.p1,
            test_args.p2)
    with raises(errors.UnitsError):
        law.calculate_magnetic_flux_density(test_args.i, test_args.n, test_args.l, nb, test_args.p2)
    with raises(errors.UnitsError):
        law.calculate_magnetic_flux_density(test_args.i, test_args.n, test_args.l, test_args.p1, nb)
