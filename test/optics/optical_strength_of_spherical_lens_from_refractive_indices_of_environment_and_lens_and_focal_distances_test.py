from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    errors,
    units,
    Quantity, prefixes
)

from symplyphysics.laws.optics import optical_strength_of_spherical_lens_from_refractive_indices_of_environment_and_lens_and_focal_distances as spherical_lens_law

# A thin, flat-convex lens is slightly immersed in water
# with its horizontal flat side (the convex surface of the lens is in the air).
# A narrow vertical beam of light falls on the lens from above,
# the axis of which passes exactly through the top of the convex surface.
# This beam is focused in water at a depth of h = 27 cm. The optical power of the lens
# in the air is D = 5 dptr. Find the refractive index of water.

# Radius of lens: R = (h * D - 1) / D = 0.07 (m)


@fixture(name="test_args")
def test_args_fixture():
    do = Quantity(1E10 * units.meters)
    di = Quantity(27 * prefixes.centi * units.meters)
    rl = Quantity(7 * prefixes.centi * units.meters)
    nm = 1
    Args = namedtuple("Args", ["do", "di", "rl", "nm", ])
    return Args(do=do, di=di, rl=rl, nm=nm)


def test_basic_refraction_index(test_args):
    result = spherical_lens_law.calculate_refraction_index_lens(test_args.do, test_args.di, test_args.rl, test_args.nm)
    assert result == approx(1.35, 0.001)


def test_bad_distance(test_args):
    db = Quantity(1 * units.coulomb)

    with raises(errors.UnitsError):
        spherical_lens_law.calculate_refraction_index_lens(db, test_args.di, test_args.rl, test_args.nm)
    with raises(TypeError):
        spherical_lens_law.calculate_refraction_index_lens(100, test_args.di, test_args.rl, test_args.nm)

    with raises(errors.UnitsError):
        spherical_lens_law.calculate_refraction_index_lens(test_args.do, db, test_args.rl, test_args.nm)
    with (raises(TypeError)):
        spherical_lens_law.calculate_refraction_index_lens(test_args.do, 100, test_args.rl, test_args.nm)


def test_bad_radius(test_args):
    rb = Quantity(1 * units.coulomb)

    with raises(errors.UnitsError):
        spherical_lens_law.calculate_refraction_index_lens(test_args.do, test_args.di, rb, test_args.nm)
    with raises(TypeError):
        spherical_lens_law.calculate_refraction_index_lens(test_args.do, test_args.di, 100, test_args.nm)
