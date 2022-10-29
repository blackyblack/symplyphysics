from collections import namedtuple
from pytest import approx, fixture

from math import pi
from sympy import sin, cos
from sympy.vector import CoordSys3D, ParametricRegion, ImplicitRegion
from symplyphysics import (
    units, convert_to, SI, expr_to_quantity, symbols
)
from symplyphysics.definitions import circulation_is_integral_along_curve as circulation_def

@fixture
def test_args():
    C = CoordSys3D('C')
    theta = symbols('theta')
    force_unit = units.Quantity('force_unit')
    SI.set_quantity_dimension(force_unit, units.force)
    SI.set_quantity_scale_factor(force_unit, 1 * units.newton)
    radius_unit = units.Quantity('radius_unit')
    SI.set_quantity_dimension(radius_unit, units.length)
    SI.set_quantity_scale_factor(radius_unit, 1 * units.meter)
    # field is a field of gravitational forces, force is directed down by the Y coordinate
    # field is (0, -1 * G * m * M / y**2)
    # G * m * M = k * force * length**2 / mass**2
    # let k = 1
    field = 0 * C.i + -1 * force_unit * radius_unit**2 / C.y**2 * C.j

    Args = namedtuple('Args', ['C', 'theta', 'force_unit', 'radius_unit', 'field'])
    return Args(C=C, theta=theta, force_unit=force_unit, radius_unit=radius_unit, field=field)


def test_basic_circulation(test_args):
    field = test_args.C.y * test_args.C.i + test_args.C.z * test_args.C.k + test_args.C.x * test_args.C.k
    curve = ParametricRegion((cos(test_args.theta), sin(test_args.theta)), (test_args.theta, 0, pi/2))
    result_expr = circulation_def.calculate_circulation(field, curve)
    assert result_expr.evalf(4) == approx(-pi/4, 0.001)

def test_two_parameters_flux(test_args):
    phi = symbols('phi')
    field = test_args.C.x**3 * test_args.C.i + test_args.C.y**3 * test_args.C.j + test_args.C.z**3 * test_args.C.k
    sphere = ParametricRegion(
        (4 * sin(phi) * cos(test_args.theta), 4 * sin(phi) * sin(test_args.theta), 4 * cos(phi)),
        (phi, 0, pi),
        (test_args.theta, 0, 2 * pi))
    # technically it is flux, not circulation
    result_expr = circulation_def.calculate_circulation(field, sphere)
    assert result_expr.evalf(4) == approx(12288 * pi/5, 0.001)

def test_two_parameters_circulation(test_args):
    x, y = symbols('x y')
    field = test_args.C.y * test_args.C.i + -1 * test_args.C.x * test_args.C.j
    ellipse = ParametricRegion((3 * cos(test_args.theta), 3 * sin(test_args.theta)), (test_args.theta, 0, 2 * pi))
    result_expr = circulation_def.calculate_circulation(field, ellipse)
    assert result_expr.evalf(4) == approx(-18 * pi, 0.001)
    #HACK: implicit regions are not really supported by sympy. Avoid using them.
    ellipse_implicit = ImplicitRegion((x, y), x**2 + y**2 - 9)
    result_expr = circulation_def.calculate_circulation(field, ellipse_implicit)
    assert result_expr.evalf(4) == approx(-18 * pi, 0.001)

def test_force_circulation(test_args):
    # trajectory is linear: y = x
    # parametrize by theta
    #HACK: gravitational force is undefined at 0 distance, use any non-zero value
    trajectory = ParametricRegion((test_args.theta, test_args.theta), (test_args.theta, 1 * test_args.radius_unit, 2 * test_args.radius_unit))

    result_expr = circulation_def.calculate_circulation(test_args.field, trajectory)
    result = expr_to_quantity(result_expr, 'force_work')
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.energy)
    result_work = convert_to(result, units.joule).subs({units.joule: 1}).evalf(2)
    assert result_work == approx(-0.5, 0.01)

def test_force_circulation_horizontal(test_args):
    # trajectory is horizontal line: y = 5
    trajectory_horizontal = ParametricRegion((test_args.theta, 5 * test_args.radius_unit), (test_args.theta, 1 * test_args.radius_unit, 2 * test_args.radius_unit))
    result_expr = circulation_def.calculate_circulation(test_args.field, trajectory_horizontal)
    assert result_expr == 0

def test_force_circulation_horizontal_up(test_args):
    # trajectory is vertical line: x = 5
    trajectory_vertical = ParametricRegion((5 * test_args.radius_unit, test_args.theta), (test_args.theta, 1 * test_args.radius_unit, 2 * test_args.radius_unit))
    result_expr = circulation_def.calculate_circulation(test_args.field, trajectory_vertical)
    result = expr_to_quantity(result_expr, 'force_work_vertical_up')
    result_work = convert_to(result, units.joule).subs({units.joule: 1}).evalf(2)
    assert result_work == approx(-0.5, 0.01)

def test_force_circulation_horizontal_down(test_args):
    # trajectory is vertical line, but with down direction: x = 6
    trajectory_vertical = ParametricRegion((6 * test_args.radius_unit, test_args.theta), (test_args.theta, 2 * test_args.radius_unit, 1 * test_args.radius_unit))
    result_expr = circulation_def.calculate_circulation(test_args.field, trajectory_vertical)
    result = expr_to_quantity(result_expr, 'force_work_vertical_down')
    result_work = convert_to(result, units.joule).subs({units.joule: 1}).evalf(2)
    assert result_work == approx(0.5, 0.01)

