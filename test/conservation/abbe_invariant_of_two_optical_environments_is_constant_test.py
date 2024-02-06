from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_approx,
    errors,
    units,
    Quantity,
    prefixes,
)
from symplyphysics.laws.conservation import abbe_invariant_of_two_optical_environments_is_constant as abbe_conservation_law

# A thin, flat-convex lens is slightly immersed in water
# with its horizontal flat side (the convex surface of the lens is in the air).
# A narrow vertical beam of light falls on the lens from above,
# the axis of which passes exactly through the top of the convex surface.
# This beam is focused in water at a depth of h = 27 cm. The optical power of the lens
# in the air is D = 5 dptr. Find the refractive index of water.

# Radius of lens: R = (h * D - 1) / D = 0.07 (m)
# Refraction index in air: n = 1

# Refraction index in lens: n' = n ( ((1/oo) - (1/R)) / ((1/h) - (1/R)) ) = 1.35


@fixture(name="test_args")
def test_args_fixture():
    r = Quantity(7 * prefixes.centi * units.meters)
    n = 1.0
    a = Quantity(1e10 * units.meters)
    b = Quantity(27 * prefixes.centi * units.meters)

    Args = namedtuple("Args", ["a", "b", "n", "r"])
    return Args(a=a, b=b, n=n, r=r)


def test_basic_conservation(test_args):
    result = abbe_conservation_law.calculate_refraction_index_lens(test_args.a, test_args.b,
        test_args.r, test_args.n)
    assert_approx(result, 1.35)


def test_bad_distance(test_args):
    db = Quantity(1 * units.coulomb)

    with raises(errors.UnitsError):
        abbe_conservation_law.calculate_refraction_index_lens(db, test_args.b, test_args.r,
            test_args.n)
    with raises(TypeError):
        abbe_conservation_law.calculate_refraction_index_lens(100, test_args.b, test_args.r,
            test_args.n)

    with raises(errors.UnitsError):
        abbe_conservation_law.calculate_refraction_index_lens(test_args.a, db, test_args.r,
            test_args.n)
    with raises(TypeError):
        abbe_conservation_law.calculate_refraction_index_lens(test_args.a, 100, test_args.r,
            test_args.n)


def test_bad_radius(test_args):
    rb = Quantity(1 * units.coulomb)

    with raises(errors.UnitsError):
        abbe_conservation_law.calculate_refraction_index_lens(
            test_args.a,
            test_args.b,
            rb,
            test_args.n,
        )
    with raises(TypeError):
        abbe_conservation_law.calculate_refraction_index_lens(
            test_args.a,
            test_args.b,
            100,
            test_args.n,
        )
