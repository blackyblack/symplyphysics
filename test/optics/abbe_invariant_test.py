from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    errors,
    units,
    convert_to,
    Quantity,
    SI, prefixes,
)
from symplyphysics.laws.optics import abbe_invariant as abbe_invariant_law

# A thin, flat-convex lens is slightly immersed in water
# with its horizontal flat side (the convex surface of the lens is in the air).
# A narrow vertical beam of light falls on the lens from above,
# the axis of which passes exactly through the top of the convex surface.
# This beam is focused in water at a depth of h = 27 cm. The optical power of the lens
# in the air is D = 5 dptr. Find the refractive index of water.

# Radius of lens: R = (h * D - 1) / D = 0.07 (m)
# Refraction index in lens: n' = h D = 1.35
# Abbe's invariant in lens: Q' = n' ((1 / h) - (1 / R)) = -14.286 (1/m)


@fixture(name="test_args")
def test_args_fixture():
    n = 1.35
    r = Quantity(7 * prefixes.centi * units.meters)
    h = Quantity(27 * prefixes.centi * units.meters)
    Args = namedtuple("Args", ["h", "n", "r"])
    return Args(h=h, n=n, r=r)


def test_basic_conservation(test_args):
    result = abbe_invariant_law.calculate_abbe_invariant(test_args.h, test_args.r, test_args.n)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, 1 / units.length)
    result_ = convert_to(result, 1 / units.meter).evalf(2)
    assert result_ == approx(-14.286, 0.01)


def test_bad_distance(test_args):
    hb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        abbe_invariant_law.calculate_abbe_invariant(hb, test_args.r, test_args.n)
    with raises(TypeError):
        abbe_invariant_law.calculate_abbe_invariant(100, test_args.r, test_args.n)
        

def test_bad_radius(test_args):
    rb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        abbe_invariant_law.calculate_abbe_invariant(test_args.h, rb, test_args.n)
    with raises(errors.UnitsError):
        abbe_invariant_law.calculate_abbe_invariant(test_args.h, 100, test_args.n)
