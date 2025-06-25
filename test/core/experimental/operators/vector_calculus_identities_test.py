from typing import Sequence, Optional

from sympy import Function as SymFunction, sqrt, S, cosh, cos, sinh, sin, exp, Expr

from symplyphysics.core.symbols.symbols import BasicSymbol

from symplyphysics import Symbol, units
from symplyphysics.core.experimental.coordinate_systems import (CartesianCoordinateSystem,
    CylindricalCoordinateSystem, SphericalCoordinateSystem, BaseCoordinateSystem, CoordinateScalar,
    CoordinateVector)
from symplyphysics.core.experimental.operators import (VectorGradient, VectorDivergence, VectorCurl,
    VectorLaplacian)
from symplyphysics.core.experimental.solvers import vector_equals


def check_div_curl(system: BaseCoordinateSystem) -> bool:
    """Identity `div(curl(v)) = 0` holds for any vector `v`."""

    a, b, c = system.base_scalars

    components = (
        SymFunction("A")(a, b, c),
        SymFunction("B")(a, b, c),
        SymFunction("C")(a, b, c),
    )

    point = BasicSymbol("P")

    vector = CoordinateVector(components, system, point)
    rot = VectorCurl(vector)
    div_rot = VectorDivergence(rot).expand().simplify()

    return div_rot == 0


def check_curl_grad(system: BaseCoordinateSystem) -> bool:
    """Identity `curl(grad(f)) = 0` hold for any scalar `f`."""

    a, b, c = system.base_scalars

    point = BasicSymbol("P")

    scalar = CoordinateScalar(SymFunction("f")(a, b, c), system, point)

    grad = VectorGradient(scalar)
    rot_grad = VectorCurl(grad)

    return rot_grad == 0


def check_laplacian_is_div_grad(system: BaseCoordinateSystem) -> bool:
    """Identity `Laplace(f) = div(grad(f))` hold for any scalar `f`."""

    a, b, c = system.base_scalars

    point = BasicSymbol("P")

    scalar = CoordinateScalar(SymFunction("f")(a, b, c), system, point)

    div_grad = VectorDivergence(VectorGradient(scalar))

    laplace = VectorLaplacian(scalar)

    return vector_equals(div_grad, laplace)


class ParabolicCoordinateSystem2D(BaseCoordinateSystem):  # pylint: disable=abstract-method

    def generate_base_scalars(self) -> tuple[Symbol, Symbol, Symbol]:
        return (
            Symbol("sigma", sqrt(units.length), display_latex="\\sigma", real=True),
            Symbol("tau", sqrt(units.length), display_latex="\\tau", real=True),
            Symbol("z", units.length, real=True),
        )

    def cartesian_transform(self, base_scalars: Optional[Sequence[Expr]] = None) -> Sequence[Expr]:
        sigma, tau, z = base_scalars or self.base_scalars

        x = sigma * tau
        y = (tau**2 - sigma**2) / 2

        return x, y, z

    def generate_lame_coefficients(self) -> tuple[Expr, Expr, Expr]:
        sigma, tau, _ = self.base_scalars

        h = sqrt(sigma**2 + tau**2)

        return h, h, S.One


class ParabolicCoordinateSystem3D(BaseCoordinateSystem):  # pylint: disable=abstract-method

    def generate_base_scalars(self) -> tuple[Symbol, Symbol, Symbol]:
        return (
            Symbol("sigma", sqrt(units.length), real=True),
            Symbol("tau", sqrt(units.length), real=True),
            Symbol("phi", units.length, display_latex="\\varphi", real=True),
        )

    def cartesian_transform(self, base_scalars: Optional[Sequence[Expr]] = None) -> Sequence[Expr]:
        sigma, tau, phi = base_scalars or self.base_scalars

        x = sigma * tau * cos(phi)
        y = sigma * tau * sin(phi)
        z = (tau**2 - sigma**2) / 2

        return x, y, z

    def generate_lame_coefficients(self) -> tuple[Expr, Expr, Expr]:
        sigma, tau, _ = self.base_scalars

        h1 = sqrt(sigma**2 + tau**2)
        h2 = sigma * tau

        return h1, h1, h2


class SixSphereCoordinateSystem(BaseCoordinateSystem):  # pylint: disable=abstract-method

    def generate_base_scalars(self) -> tuple[Symbol, Symbol, Symbol]:
        return (
            Symbol("u", 1 / units.length, real=True),
            Symbol("v", 1 / units.length, real=True),
            Symbol("w", 1 / units.length, real=True),
        )

    def cartesian_transform(self, base_scalars: Optional[Sequence[Expr]] = None) -> Sequence[Expr]:
        u, v, w = base_scalars or self.base_scalars

        s = u**2 + v**2 + w**2

        return u / s, v / s, w / s

    def generate_lame_coefficients(self) -> tuple[Expr, Expr, Expr]:
        u, v, w = self.base_scalars

        s = u**2 + v**2 + w**2

        return 1 / s, 1 / s, 1 / s


class EllipticalCylindricalCoordinateSystem(BaseCoordinateSystem):  # pylint: disable=abstract-method
    _factor = Symbol("a", units.length, real=True)

    @property
    def factor(self) -> Symbol:
        return self._factor

    def generate_base_scalars(self) -> tuple[Symbol, Symbol, Symbol]:
        return (
            Symbol("mu", nonnegative=True, display_latex="\\mu"),
            Symbol("nu", real=True, display_latex="\\nu"),
            Symbol("z", real=True),
        )

    def cartesian_transform(self, base_scalars: Optional[Sequence[Expr]] = None) -> Sequence[Expr]:
        mu, nu, z = base_scalars or self.base_scalars

        x = self.factor * cosh(mu) * cos(nu)
        y = self.factor * sinh(mu) * sin(nu)

        return x, y, z

    def generate_lame_coefficients(self) -> tuple[Expr, Expr, Expr]:
        mu, nu, _ = self.base_scalars

        h = self.factor * sqrt(sinh(mu)**2 + sin(nu)**2)

        return h, h, S.One


class LogPolarCoordinateSystem(BaseCoordinateSystem):  # pylint: disable=abstract-method
    _factor = Symbol("a", units.length, real=True)

    @property
    def factor(self) -> Symbol:
        return self._factor

    def generate_base_scalars(self) -> tuple[Symbol, Symbol, Symbol]:
        return (
            Symbol("rho", real=True),
            Symbol("phi", real=True),
            Symbol("z", real=True),
        )

    def cartesian_transform(self, base_scalars: Optional[Sequence[Expr]] = None) -> Sequence[Expr]:
        rho, phi, z = base_scalars or self.base_scalars

        x = self.factor * exp(rho) * cos(phi)
        y = self.factor * exp(rho) * sin(phi)

        return x, y, z

    def generate_lame_coefficients(self) -> tuple[Expr, Expr, Expr]:
        rho, _, _ = self.base_scalars

        h = abs(self.factor) * exp(rho)

        return h, h, S.One


def test_cartesian_system() -> None:
    cartesian = CartesianCoordinateSystem()

    assert check_div_curl(cartesian)
    assert check_curl_grad(cartesian)
    assert check_laplacian_is_div_grad(cartesian)


def test_cylindrical_system() -> None:
    cylindrical = CylindricalCoordinateSystem()

    assert check_div_curl(cylindrical)
    assert check_curl_grad(cylindrical)
    assert check_laplacian_is_div_grad(cylindrical)


def test_spherical_system() -> None:
    spherical = SphericalCoordinateSystem()

    assert check_div_curl(spherical)
    assert check_curl_grad(spherical)
    assert check_laplacian_is_div_grad(spherical)


def test_parabolic_2d() -> None:
    parabolic3d = ParabolicCoordinateSystem2D()
    assert check_div_curl(parabolic3d)
    assert check_curl_grad(parabolic3d)
    assert check_laplacian_is_div_grad(parabolic3d)


def test_parabolic_3d() -> None:
    parabolic3d = ParabolicCoordinateSystem3D()
    assert check_div_curl(parabolic3d)
    assert check_curl_grad(parabolic3d)
    assert check_laplacian_is_div_grad(parabolic3d)


def test_six_sphere() -> None:
    sphere6 = SixSphereCoordinateSystem()
    assert check_div_curl(sphere6)
    assert check_curl_grad(sphere6)
    assert check_laplacian_is_div_grad(sphere6)


def test_elliptical_cylindrical() -> None:
    elcyl = EllipticalCylindricalCoordinateSystem()
    assert check_div_curl(elcyl)
    assert check_curl_grad(elcyl)
    assert check_laplacian_is_div_grad(elcyl)


def test_logpolar() -> None:
    logpolar = LogPolarCoordinateSystem()
    assert check_div_curl(logpolar)
    assert check_curl_grad(logpolar)
    assert check_laplacian_is_div_grad(logpolar)
