from collections import namedtuple
from pytest import approx, fixture, raises
from sympy import oo

from symplyphysics import (
    errors,
    units,
    Quantity, prefixes
)

from symplyphysics.laws.optics import the_formula_of_the_spherical_refractive_surface as spherical_lens_law

# A thin, flat-convex lens is slightly immersed in water
# with its horizontal flat side (the convex surface of the lens is in the air).
# A narrow vertical beam of light falls on the lens from above,
# the axis of which passes exactly through the top of the convex surface.
# This beam is focused in water at a depth of h = 27 cm. The optical power of the lens
# in the air is D = 5 dptr. Find the refractive index of water.

# Radius of lens: R = (h * D - 1) / D = 0.07 (m)


@fixture(name="test_args")
def test_args_fixture():
    do = Quantity(oo * units.meters)
    di = Quantity(27 * prefixes.centi * units.meters)
    nm = 1
    rl = Quantity(7 * prefixes.centi * units.meter)
    Args = namedtuple("Args", ["do", "di", "nm", "rl"])
    return Args(do=do, di=di, nm=nm, rl=rl)


def test_basic_refraction_index(test_args):
    result = spherical_lens_law.calculate_refraction_index_lens(test_args.do, test_args.di, test_args.nm, test_args.rl)
    assert result == approx(1.35, 0.001)


def test_bad_distance(test_args):
    db = Quantity(1 * units.coulomb)

    with raises(errors.UnitsError):
        spherical_lens_law.calculate_refraction_index_lens(db, test_args.di, test_args.nm, test_args.rl)
    with raises(TypeError):
        spherical_lens_law.calculate_refraction_index_lens(100, test_args.di, test_args.nm, test_args.rl)

    with raises(errors.UnitsError):
        spherical_lens_law.calculate_refraction_index_lens(test_args.do, db, test_args.nm, test_args.rl)
    with (raises(TypeError)):
        spherical_lens_law.calculate_refraction_index_lens(test_args.do, 100, test_args.nm, test_args.rl)


def test_bad_radius(test_args):
    rb = Quantity(1 * units.coulomb)

    with raises(errors.UnitsError):
        spherical_lens_law.calculate_refraction_index_lens(test_args.do, test_args.di, test_args.nm, rb)
    with raises(TypeError):
        spherical_lens_law.calculate_refraction_index_lens(test_args.do, test_args.di, test_args.nm, 100)
