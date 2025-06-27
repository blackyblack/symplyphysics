from collections import namedtuple
from pytest import fixture
from sympy import cos, pi, sin, Symbol as SymSymbol
from symplyphysics import assert_equal, units
from symplyphysics.laws.dynamics.vector import mechanical_work_is_line_integral_of_force as law
from symplyphysics.core.expr_comparisons import expr_equals

from symplyphysics.core.experimental.coordinate_systems import (CoordinateVector, CARTESIAN,
    AppliedPoint)
from symplyphysics.core.experimental.coordinate_systems.curve import Curve

Args = namedtuple("Args", "f0 fxy c")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    x, y, _ = CARTESIAN.base_scalars

    f0 = CoordinateVector([units.newton, 0, 0], CARTESIAN)
    fxy = CoordinateVector([
        units.newton * sin(y / units.meter),
        0,
        units.newton * cos(x / units.meter),
    ], CARTESIAN)

    t = SymSymbol("t", real=True)
    parametrization = AppliedPoint([
        -1 * units.meter * cos(t),
        2 * units.meter * sin(t),
        0,
    ], CARTESIAN)
    c = Curve(t, parametrization)

    return Args(f0=f0, fxy=fxy, c=c)


def test_law(test_args: Args) -> None:
    result_0 = law.calculate_work(test_args.f0, test_args.c, 0, pi)
    assert_equal(result_0, 2 * units.joule)

    result_xy = law.calculate_work(test_args.fxy, test_args.c, 0, pi)
    assert_equal(result_xy, 1.81 * units.joule, relative_tolerance=2e-2)


def test_law_closed_curve(test_args: Args) -> None:
    result_0 = law.calculate_work(test_args.f0, test_args.c, 0, 2 * pi)
    assert_equal(result_0, 0)  # the force is conservative

    result_xy = law.calculate_work(test_args.fxy, test_args.c, 0, 2 * pi)
    # the force isn't conservative; also note that it's twice the value in the test above
    assert_equal(result_xy, 3.62 * units.joule, relative_tolerance=2e-3)


def test_orthogonal_force(test_args: Args) -> None:
    t = SymSymbol("t")
    parametrization = AppliedPoint([0, units.meter * t, 0], CARTESIAN)
    c = Curve(t, parametrization)

    result_0 = law.law.rhs.subs({
        law.force(law.position_vector): test_args.f0,
        law.curve: c,
    }).doit().doit()
    assert expr_equals(result_0, 0)

    result_xy = law.law.rhs.subs({
        law.force(law.position_vector): test_args.fxy,
        law.curve: c,
    }).doit().doit()
    assert expr_equals(result_xy, 0)
