from dataclasses import dataclass
from pytest import raises, fixture
from sympy import Symbol as SymSymbol, cos, sin
from sympy.physics.units.systems.si import dimsys_SI
from symplyphysics import symbols, units, dimensionless, Quantity
from symplyphysics.core.errors import UnitsError
from symplyphysics.core.dimensions.collect_expression import collect_expression_and_dimension

from symplyphysics.core.experimental.vectors import (VectorSymbol, VectorFunction, VectorDot,
    VectorCross, VectorNorm, VectorMixedProduct, vector_diff)
from symplyphysics.core.experimental.coordinate_systems import (CoordinateScalar, CoordinateVector,
    QuantityCoordinateVector, CartesianCoordinateSystem, CylindricalCoordinateSystem,
    SphericalCoordinateSystem)
from symplyphysics.core.experimental.operators import (VectorGradient, VectorDivergence, VectorCurl,
    VectorLaplacian)


@dataclass(frozen=True, kw_only=True)
class Args:
    cart: CartesianCoordinateSystem
    cyl: CylindricalCoordinateSystem
    sph: SphericalCoordinateSystem
    p: SymSymbol


@fixture(name="test_args")
def test_args_fixture() -> Args:
    cart = CartesianCoordinateSystem()
    cyl = CylindricalCoordinateSystem()
    sph = SphericalCoordinateSystem()
    p = SymSymbol("P")

    return Args(cart=cart, cyl=cyl, sph=sph, p=p)


def test_vector_symbol() -> None:
    expr = VectorSymbol("v")
    dim = collect_expression_and_dimension(expr)[1]
    assert dimsys_SI.is_dimensionless(dim)

    expr = VectorSymbol("v", units.length)
    dim = collect_expression_and_dimension(expr)[1]
    assert dimsys_SI.equivalent_dims(dim, expr.dimension)

    expr = VectorSymbol("v", units.momentum / units.time)
    dim = collect_expression_and_dimension(expr)[1]
    assert dimsys_SI.equivalent_dims(dim, expr.dimension)


def test_vector_function() -> None:
    func = VectorFunction("f")
    dim = collect_expression_and_dimension(func)[1]
    assert dimsys_SI.is_dimensionless(dim)

    func = VectorFunction("f", dimension=units.length)
    dim = collect_expression_and_dimension(func)[1]
    assert dimsys_SI.equivalent_dims(dim, func.dimension)

    func = VectorFunction("f", dimension=units.momentum / units.time)
    dim = collect_expression_and_dimension(func)[1]
    assert dimsys_SI.equivalent_dims(dim, func.dimension)


def test_applied_vector_function() -> None:
    t = symbols.time

    expr = VectorFunction("f")(t)
    dim = collect_expression_and_dimension(expr)[1]
    assert dimsys_SI.is_dimensionless(dim)

    expr = VectorFunction("f", dimension=units.length)(t)
    dim = collect_expression_and_dimension(expr)[1]
    assert dimsys_SI.equivalent_dims(dim, units.length)

    expr = VectorFunction("f", dimension=units.charge)(t)
    dim = collect_expression_and_dimension(expr)[1]
    assert dimsys_SI.equivalent_dims(dim, units.charge)


def test_vector_mul() -> None:
    v = VectorSymbol("v")

    expr = -1 * v
    dim = collect_expression_and_dimension(expr)[1]
    assert dimsys_SI.is_dimensionless(dim)

    f = VectorSymbol("f", units.force)
    expr = 4 * f
    dim = collect_expression_and_dimension(expr)[1]
    assert dimsys_SI.equivalent_dims(dim, units.force)

    expr = VectorNorm(f) * f
    dim = collect_expression_and_dimension(expr)[1]
    assert dimsys_SI.equivalent_dims(dim, units.force**2)

    expr = 4 * units.second * symbols.position * f
    dim = collect_expression_and_dimension(expr)[1]
    assert dimsys_SI.equivalent_dims(dim, units.force * units.time * units.length)

    expr = (VectorNorm(f) + 4 * units.newton) * v
    dim = collect_expression_and_dimension(expr)[1]
    assert dimsys_SI.equivalent_dims(dim, units.force)


def test_vector_add() -> None:
    v = VectorSymbol("v")
    w = VectorSymbol("w")

    expr = v + 0
    dim = collect_expression_and_dimension(expr)[1]
    assert dimsys_SI.is_dimensionless(dim)

    expr = 0 + w
    dim = collect_expression_and_dimension(expr)[1]
    assert dimsys_SI.is_dimensionless(dim)

    expr = v + w
    dim = collect_expression_and_dimension(expr)[1]
    assert dimsys_SI.is_dimensionless(dim)

    f = VectorSymbol("f", units.force)
    g = VectorSymbol("g", units.force)
    h = VectorSymbol("h", units.force)

    expr = 0 + g
    dim = collect_expression_and_dimension(expr)[1]
    assert dimsys_SI.equivalent_dims(dim, units.force)

    expr = f + 0
    dim = collect_expression_and_dimension(expr)[1]
    assert dimsys_SI.equivalent_dims(dim, units.force)

    expr = f + g
    dim = collect_expression_and_dimension(expr)[1]
    assert dimsys_SI.equivalent_dims(dim, units.force)

    expr = f + g + h
    dim = collect_expression_and_dimension(expr)[1]
    assert dimsys_SI.equivalent_dims(dim, units.force)

    expr = f + g + 0
    dim = collect_expression_and_dimension(expr)[1]
    assert dimsys_SI.equivalent_dims(dim, units.force)

    expr = VectorCross(v, f) + VectorCross(g, w)
    dim = collect_expression_and_dimension(expr)[1]
    assert dimsys_SI.equivalent_dims(dim, units.force)

    expr = 1 + f
    with raises(UnitsError):
        collect_expression_and_dimension(expr)

    expr = v + f
    with raises(UnitsError):
        collect_expression_and_dimension(expr)

    expr = v + w + f + g
    with raises(UnitsError):
        collect_expression_and_dimension(expr)


def test_vector_dot() -> None:
    v = VectorSymbol("v", units.velocity)
    p = VectorSymbol("p", units.momentum)
    r = VectorSymbol("r", units.length)

    expr = VectorDot(v, p)
    dim = collect_expression_and_dimension(expr)[1]
    assert dimsys_SI.equivalent_dims(dim, units.velocity * units.momentum)

    expr = VectorDot(p, v)
    dim = collect_expression_and_dimension(expr)[1]
    assert dimsys_SI.equivalent_dims(dim, units.velocity * units.momentum)

    expr = VectorDot(0, p, evaluate=False)
    dim = collect_expression_and_dimension(expr)[1]
    assert dimsys_SI.equivalent_dims(dim, dimensionless)

    expr = VectorDot(p, 0, evaluate=False)
    dim = collect_expression_and_dimension(expr)[1]
    assert dimsys_SI.equivalent_dims(dim, dimensionless)

    expr = VectorDot(VectorCross(v, p), r)
    dim = collect_expression_and_dimension(expr)[1]
    assert dimsys_SI.equivalent_dims(dim, units.velocity * units.momentum * units.length)


def test_vector_norm() -> None:
    expr = VectorNorm(0, evaluate=False)
    dim = collect_expression_and_dimension(expr)[1]
    assert dimsys_SI.is_dimensionless(dim)

    v = VectorSymbol("v", units.velocity)
    w = VectorSymbol("w", units.velocity)

    expr = VectorNorm(v)
    dim = collect_expression_and_dimension(expr)[1]
    assert dimsys_SI.equivalent_dims(dim, units.velocity)

    expr = VectorNorm(v + 0)
    dim = collect_expression_and_dimension(expr)[1]
    assert dimsys_SI.equivalent_dims(dim, units.velocity)

    expr = VectorNorm(0 + v)
    dim = collect_expression_and_dimension(expr)[1]
    assert dimsys_SI.equivalent_dims(dim, units.velocity)

    expr = VectorNorm(v + v)
    dim = collect_expression_and_dimension(expr)[1]
    assert dimsys_SI.equivalent_dims(dim, units.velocity)

    expr = VectorNorm(v + w)
    dim = collect_expression_and_dimension(expr)[1]
    assert dimsys_SI.equivalent_dims(dim, units.velocity)


def test_vector_mixed_product() -> None:
    a = VectorSymbol("a", units.length)
    b = VectorSymbol("b", units.force)
    c = VectorSymbol("c", units.momentum)
    d = VectorSymbol("d", units.momentum)

    expr = VectorMixedProduct(0, b, c, evaluate=False)
    dim = collect_expression_and_dimension(expr)[1]
    assert dimsys_SI.equivalent_dims(dim, dimensionless)

    expr = VectorMixedProduct(a, 0, c, evaluate=False)
    dim = collect_expression_and_dimension(expr)[1]
    assert dimsys_SI.equivalent_dims(dim, dimensionless)

    expr = VectorMixedProduct(a, b, 0, evaluate=False)
    dim = collect_expression_and_dimension(expr)[1]
    assert dimsys_SI.equivalent_dims(dim, dimensionless)

    expr = VectorMixedProduct(a, b, c)
    dim = collect_expression_and_dimension(expr)[1]
    assert dimsys_SI.equivalent_dims(dim, units.length * units.force * units.momentum)

    expr = VectorMixedProduct(a, b, c + d)
    dim = collect_expression_and_dimension(expr)[1]
    assert dimsys_SI.equivalent_dims(dim, units.length * units.force * units.momentum)

    expr = VectorMixedProduct(b, c, a)
    dim = collect_expression_and_dimension(expr)[1]
    assert dimsys_SI.equivalent_dims(dim, units.length * units.force * units.momentum)

    expr = VectorMixedProduct(b, c + d, a)
    dim = collect_expression_and_dimension(expr)[1]
    assert dimsys_SI.equivalent_dims(dim, units.length * units.force * units.momentum)

    expr = VectorMixedProduct(c, a, b)
    dim = collect_expression_and_dimension(expr)[1]
    assert dimsys_SI.equivalent_dims(dim, units.length * units.force * units.momentum)

    expr = VectorMixedProduct(c + d, a, b)
    dim = collect_expression_and_dimension(expr)[1]
    assert dimsys_SI.equivalent_dims(dim, units.length * units.force * units.momentum)


def test_vector_cross() -> None:
    v = VectorSymbol("v", units.velocity)
    w = VectorSymbol("w", units.velocity)
    p = VectorSymbol("p", units.momentum)
    r = VectorSymbol("r", units.length)

    expr = VectorCross(0, p, evaluate=False)
    dim = collect_expression_and_dimension(expr)[1]
    assert dimsys_SI.equivalent_dims(dim, dimensionless)

    expr = VectorCross(v, 0, evaluate=False)
    dim = collect_expression_and_dimension(expr)[1]
    assert dimsys_SI.equivalent_dims(dim, dimensionless)

    expr = VectorCross(v, p)
    dim = collect_expression_and_dimension(expr)[1]
    assert dimsys_SI.equivalent_dims(dim, units.velocity * units.momentum)

    expr = VectorCross(p, v)
    dim = collect_expression_and_dimension(expr)[1]
    assert dimsys_SI.equivalent_dims(dim, units.velocity * units.momentum)

    expr = VectorCross(p, v + w)
    dim = collect_expression_and_dimension(expr)[1]
    assert dimsys_SI.equivalent_dims(dim, units.velocity * units.momentum)

    expr = VectorCross(v + w, p)
    dim = collect_expression_and_dimension(expr)[1]
    assert dimsys_SI.equivalent_dims(dim, units.velocity * units.momentum)

    expr = VectorCross(VectorCross(v, p), r)
    dim = collect_expression_and_dimension(expr)[1]
    assert dimsys_SI.equivalent_dims(dim, units.velocity * units.momentum * units.length)

    expr = VectorCross(r, VectorCross(v, p))
    dim = collect_expression_and_dimension(expr)[1]
    assert dimsys_SI.equivalent_dims(dim, units.velocity * units.momentum * units.length)


def test_vector_derivative() -> None:
    x = symbols.position
    t = symbols.time
    f = VectorFunction("F", (x, t), dimension=units.force)

    expr = vector_diff(f(x, t), t)
    dim = collect_expression_and_dimension(expr)[1]
    assert dimsys_SI.equivalent_dims(dim, units.force / units.time)

    expr = vector_diff(f(x, t), x)
    dim = collect_expression_and_dimension(expr)[1]
    assert dimsys_SI.equivalent_dims(dim, units.force / units.length)

    expr = vector_diff(f(x, t), t, x)
    dim = collect_expression_and_dimension(expr)[1]
    assert dimsys_SI.equivalent_dims(dim, units.force / units.time / units.length)

    expr = vector_diff(f(x, t), x, t)
    dim = collect_expression_and_dimension(expr)[1]
    assert dimsys_SI.equivalent_dims(dim, units.force / units.time / units.length)


def test_coordinate_scalar(test_args: Args) -> None:
    expr = CoordinateScalar(0, test_args.cart)
    dim = collect_expression_and_dimension(expr)[1]
    assert dimsys_SI.equivalent_dims(dim, dimensionless)

    expr = CoordinateScalar(1, test_args.cart)
    dim = collect_expression_and_dimension(expr)[1]
    assert dimsys_SI.equivalent_dims(dim, dimensionless)

    expr = CoordinateScalar(4 * units.meter, test_args.cart)
    dim = collect_expression_and_dimension(expr)[1]
    assert dimsys_SI.equivalent_dims(dim, units.length)


def test_coordinate_vector(test_args: Args) -> None:
    expr = CoordinateVector([0, 0, 0], test_args.cart)
    dim = collect_expression_and_dimension(expr)[1]
    assert dimsys_SI.is_dimensionless(dim)

    expr = CoordinateVector([
        0,
        Quantity(1 * units.newton),
        Quantity(1 * units.newton),
    ], test_args.cart)
    dim = collect_expression_and_dimension(expr)[1]
    assert dimsys_SI.equivalent_dims(dim, units.force)

    expr = CoordinateVector([
        0,
        0,
        Quantity(1 * units.newton),
    ], test_args.cart)
    dim = collect_expression_and_dimension(expr)[1]
    assert dimsys_SI.equivalent_dims(dim, units.force)

    components = [
        Quantity(1 * units.newton),
        Quantity(1 * units.newton),
        Quantity(1 * units.newton),
    ]

    expr = CoordinateVector(components, test_args.cart)
    dim = collect_expression_and_dimension(expr)[1]
    assert dimsys_SI.equivalent_dims(dim, units.force)

    expr = CoordinateVector(components, test_args.cyl, test_args.p)
    dim = collect_expression_and_dimension(expr)[1]
    assert dimsys_SI.equivalent_dims(dim, units.force)

    expr = CoordinateVector(components, test_args.sph, test_args.p)
    dim = collect_expression_and_dimension(expr)[1]
    assert dimsys_SI.equivalent_dims(dim, units.force)

    expr = CoordinateVector([
        Quantity(1 * units.newton),
        Quantity(1 * units.meter),
        Quantity(1 * units.newton),
    ], test_args.cart)
    with raises(UnitsError):
        collect_expression_and_dimension(expr)

    expr = CoordinateVector([
        Quantity(1 * units.newton),
        Quantity(1 * units.meter),
        Quantity(1 * units.meter),
    ], test_args.cart)
    with raises(UnitsError):
        collect_expression_and_dimension(expr)


def test_quantity_coordinate_vector(test_args: Args) -> None:
    expr = QuantityCoordinateVector([0, 0, 0], test_args.cart)
    dim = collect_expression_and_dimension(expr)[1]
    assert dimsys_SI.is_dimensionless(dim)

    expr = QuantityCoordinateVector([
        0,
        Quantity(1 * units.newton),
        Quantity(1 * units.newton),
    ], test_args.cart)
    dim = collect_expression_and_dimension(expr)[1]
    assert dimsys_SI.equivalent_dims(dim, units.force)

    expr = QuantityCoordinateVector([
        0,
        0,
        Quantity(1 * units.newton),
    ], test_args.cart)
    dim = collect_expression_and_dimension(expr)[1]
    assert dimsys_SI.equivalent_dims(dim, units.force)

    components = [
        Quantity(1 * units.newton),
        Quantity(1 * units.newton),
        Quantity(1 * units.newton),
    ]

    expr = QuantityCoordinateVector(components, test_args.cart)
    dim = collect_expression_and_dimension(expr)[1]
    assert dimsys_SI.equivalent_dims(dim, units.force)

    expr = QuantityCoordinateVector(components, test_args.cyl, test_args.p)
    dim = collect_expression_and_dimension(expr)[1]
    assert dimsys_SI.equivalent_dims(dim, units.force)

    expr = QuantityCoordinateVector(components, test_args.sph, test_args.p)
    dim = collect_expression_and_dimension(expr)[1]
    assert dimsys_SI.equivalent_dims(dim, units.force)

    with raises(UnitsError):
        _ = QuantityCoordinateVector([
            Quantity(1 * units.newton),
            Quantity(1 * units.meter),
            Quantity(1 * units.newton),
        ], test_args.cart)

    with raises(UnitsError):
        _ = QuantityCoordinateVector([
            Quantity(1 * units.newton),
            Quantity(1 * units.meter),
            Quantity(1 * units.meter),
        ], test_args.cart)


def test_vector_gradient(test_args: Args) -> None:
    x = SymSymbol("x")
    expr = VectorGradient(x, evaluate=False)
    dim = collect_expression_and_dimension(expr)[1]
    assert dimsys_SI.is_dimensionless(dim)

    x, y, z = test_args.cart.base_scalars
    s = CoordinateScalar(x**2 + y**2 + z**2, test_args.cart)
    expr = VectorGradient(s, evaluate=False)
    dim = collect_expression_and_dimension(expr)[1]
    assert dimsys_SI.equivalent_dims(dim, units.length)

    s = CoordinateScalar(x + y * z, test_args.cart)
    expr = VectorGradient(s, evaluate=False)
    with raises(UnitsError):
        collect_expression_and_dimension(expr)

    rho, phi, z = test_args.cyl.base_scalars
    s = CoordinateScalar(rho**2 + z**2 * cos(phi), test_args.cyl, test_args.p)
    expr = VectorGradient(s, evaluate=False)
    dim = collect_expression_and_dimension(expr)[1]
    assert dimsys_SI.equivalent_dims(dim, units.length)

    s = CoordinateScalar(rho + 1 / z, test_args.cyl, test_args.p)
    expr = VectorGradient(s, evaluate=False)
    with raises(UnitsError):
        collect_expression_and_dimension(expr)

    r, theta, phi = test_args.sph.base_scalars
    s = CoordinateScalar(r**2 * sin(theta + phi) * cos(theta - phi), test_args.sph, test_args.p)
    expr = VectorGradient(s, evaluate=False)
    dim = collect_expression_and_dimension(expr)[1]
    assert dimsys_SI.equivalent_dims(dim, units.length)

    s = CoordinateScalar(r**2 + sin(theta) + cos(phi), test_args.sph, test_args.p)
    expr = VectorGradient(s, evaluate=False)
    with raises(UnitsError):
        collect_expression_and_dimension(expr)


def test_vector_divergence(test_args: Args) -> None:
    f = VectorSymbol("f", units.length)
    expr = VectorDivergence(f, evaluate=False)
    dim = collect_expression_and_dimension(expr)[1]
    assert dimsys_SI.is_dimensionless(dim)

    x, y, z = test_args.cart.base_scalars
    v = CoordinateVector([x * y, y * z, z * x], test_args.cart)
    expr = VectorDivergence(v, evaluate=False)
    dim = collect_expression_and_dimension(expr)[1]
    assert dimsys_SI.equivalent_dims(dim, units.length)

    v = CoordinateVector([x * y, x, y], test_args.cart)
    expr = VectorDivergence(v, evaluate=False)
    dim = collect_expression_and_dimension(expr)[1]
    assert dimsys_SI.equivalent_dims(dim, units.length)

    v = CoordinateVector([x * y, x, z], test_args.cart)
    expr = VectorDivergence(v, evaluate=False)
    with raises(UnitsError):
        collect_expression_and_dimension(expr)

    rho, phi, z = test_args.cyl.base_scalars
    v = CoordinateVector([rho**2, rho * z * cos(phi), z**2], test_args.cyl, test_args.p)
    expr = VectorDivergence(v, evaluate=False)
    dim = collect_expression_and_dimension(expr)[1]
    assert dimsys_SI.equivalent_dims(dim, units.length)

    v = CoordinateVector([rho**3, rho * z, z**2], test_args.cyl, test_args.p)
    expr = VectorDivergence(v, evaluate=False)
    with raises(UnitsError):
        collect_expression_and_dimension(expr)

    r, theta, phi = test_args.sph.base_scalars
    v = CoordinateVector([r**2, r**2 * sin(theta), 0], test_args.sph, test_args.p)
    expr = VectorDivergence(v, evaluate=False)
    dim = collect_expression_and_dimension(expr)[1]
    assert dimsys_SI.equivalent_dims(dim, units.length)

    v = CoordinateVector([r**3, r**2 * sin(theta), 0], test_args.sph, test_args.p)
    expr = VectorDivergence(v, evaluate=False)
    with raises(UnitsError):
        collect_expression_and_dimension(expr)


def test_vector_curl(test_args: Args) -> None:
    f = VectorSymbol("f", units.length)
    expr = VectorCurl(f, evaluate=False)
    dim = collect_expression_and_dimension(expr)[1]
    assert dimsys_SI.is_dimensionless(dim)

    x, y, z = test_args.cart.base_scalars
    v = CoordinateVector([x**2 + y**2, y**2, z**3 / x], test_args.cart)
    expr = VectorCurl(v, evaluate=False)
    dim = collect_expression_and_dimension(expr)[1]
    assert dimsys_SI.equivalent_dims(dim, units.length)

    v = QuantityCoordinateVector([
        Quantity(1 * units.newton),
        Quantity(1 * units.newton),
        Quantity(1 * units.newton),
    ], test_args.cart)
    expr = VectorCurl(v, evaluate=False)
    dim = collect_expression_and_dimension(expr)[1]
    assert dimsys_SI.is_dimensionless(dim)

    v = CoordinateVector([x, y, z], test_args.cart)
    expr = VectorCurl(v, evaluate=False)
    dim = collect_expression_and_dimension(expr)[1]
    assert dimsys_SI.is_dimensionless(dim)

    rho, phi, z = test_args.cyl.base_scalars
    v = CoordinateVector([rho**2, rho * z * cos(phi), z**2], test_args.cyl, test_args.p)
    expr = VectorCurl(v, evaluate=False)
    dim = collect_expression_and_dimension(expr)[1]
    assert dimsys_SI.equivalent_dims(dim, units.length)

    v = CoordinateVector([z**2, rho * z * cos(phi), rho], test_args.cyl, test_args.p)
    expr = VectorCurl(v, evaluate=False)
    with raises(UnitsError):
        collect_expression_and_dimension(expr)

    r, theta, phi = test_args.sph.base_scalars
    v = CoordinateVector([r * sin(theta), r * cos(phi), r * sin(theta)], test_args.sph, test_args.p)
    expr = VectorCurl(v, evaluate=False)
    dim = collect_expression_and_dimension(expr)[1]
    assert dimsys_SI.is_dimensionless(dim)

    v = CoordinateVector([sin(theta), r * cos(phi), r * sin(theta)], test_args.sph, test_args.p)
    expr = VectorCurl(v, evaluate=False)
    with raises(UnitsError):
        collect_expression_and_dimension(expr)


def test_vector_laplacian(test_args: Args) -> None:
    s = SymSymbol("x")
    expr = VectorLaplacian(s, evaluate=False)
    dim = collect_expression_and_dimension(expr)[1]
    assert dimsys_SI.is_dimensionless(dim)

    x, y, z = test_args.cart.base_scalars
    s = CoordinateScalar(x**3 + y**3 + z**3, test_args.cart)
    expr = VectorLaplacian(s, evaluate=False)
    dim = collect_expression_and_dimension(expr)[1]
    assert dimsys_SI.equivalent_dims(dim, units.length)

    s = CoordinateScalar(x + y + z, test_args.cart)
    expr = VectorLaplacian(s, evaluate=False)
    dim = collect_expression_and_dimension(expr)[1]
    assert dimsys_SI.is_dimensionless(dim)

    s = CoordinateScalar(x + y**2 + z, test_args.cart)
    expr = VectorLaplacian(s, evaluate=False)
    dim = collect_expression_and_dimension(expr)[1]
    assert dimsys_SI.is_dimensionless(dim)

    s = CoordinateScalar(x + y**2 + z**3, test_args.cart)
    expr = VectorLaplacian(s, evaluate=False)
    with raises(UnitsError):
        collect_expression_and_dimension(expr)

    rho, _, z = test_args.cyl.base_scalars
    s = CoordinateScalar(sin(rho / z), test_args.cyl, test_args.p)
    expr = VectorLaplacian(s, evaluate=False)
    dim = collect_expression_and_dimension(expr)[1]
    assert dimsys_SI.equivalent_dims(dim, 1 / units.length**2)

    s = CoordinateScalar(sin(rho / z) + rho**3, test_args.cyl, test_args.p)
    expr = VectorLaplacian(s, evaluate=False)
    with raises(UnitsError):
        collect_expression_and_dimension(expr)

    r, theta, phi = test_args.sph.base_scalars
    s = CoordinateScalar(sin(theta + phi), test_args.sph, test_args.p)
    expr = VectorLaplacian(s, evaluate=False)
    dim = collect_expression_and_dimension(expr)[1]
    assert dimsys_SI.equivalent_dims(dim, 1 / units.length**2)

    s = CoordinateScalar(sin(theta) + r, test_args.sph, test_args.p)
    expr = VectorLaplacian(s, evaluate=False)
    with raises(UnitsError):
        collect_expression_and_dimension(expr)
