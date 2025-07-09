from __future__ import annotations

from typing import Optional, Any, Sequence

from sympy import Expr
from sympy.physics import units

from symplyphysics.core.symbols.symbols import BasicSymbol, Function

from ..miscellaneous import evaluate_or_global_fallback

from ..vectors import VectorSymbol, vector_diff, VectorCross
from ..coordinate_systems import CoordinateVector
from ..coordinate_systems.vector import construct_position_vector
from ..coordinate_systems.surface import Surface

INFINITESIMAL_VECTOR_AREA = VectorSymbol("dS", units.area, display_latex="d \\vec S")


class SurfaceIntegral(Expr):  # pylint: disable=too-few-public-methods
    """
    A **surface integral** is an integral where the integrand is evaluated over a surface.

    Surface integral formulas can be found
    `here <https://en.wikipedia.org/wiki/Curvilinear_coordinates#Integration_2>`__.

    ..
        NOTE: we do not distinguish yet between integrals over open and closed surfaces
    """

    def __new__(
        cls,
        expr: Expr,
        surface: Surface | BasicSymbol,
        bounds: Optional[Sequence[tuple[Any, Any]]] = None,
        *,
        evaluate: Optional[bool] = None,
    ) -> Expr:
        evaluate = evaluate_or_global_fallback(evaluate)

        if not evaluate or isinstance(surface, BasicSymbol):
            if not bounds:
                return super().__new__(cls, expr, surface)  # pylint: disable=too-many-function-args

            return super().__new__(cls, expr, surface, bounds)  # pylint: disable=too-many-function-args

        if bounds is None:
            raise ValueError("Define the surface parameters' bounds for the surface integral")

        (t10, t11), (t20, t21) = bounds

        system = surface.parametrization.system
        position_vector_components = surface.parametrization.coordinates.values()

        expr = expr.subs(system.base_scalar_subs(position_vector_components))

        r = construct_position_vector(surface.parametrization)

        t1, t2 = surface.parameters

        r0, r1, r2 = r.components

        # NOTE: a temporary fix due to VectorDerivative returning 0 when the system is
        # non-Cartesian and the differentiation symbols don't appear in the vector components
        g = Function("g")
        r = CoordinateVector([r0 + g(t1, t2), r1, r2], r.system, r.point)

        dr_dt1 = vector_diff(r, t1).subs(g(t1, t2), 0).doit()
        dr_dt2 = vector_diff(r, t2).subs(g(t1, t2), 0).doit()

        n = VectorCross(dr_dt1, dr_dt2).simplify()

        expr = expr.subs(INFINITESIMAL_VECTOR_AREA, n)

        for f, q in zip(system.base_scalar_functions, position_vector_components):
            expr = expr.subs({f(t1): q, f(t2): q})

        return surface.to_integral(expr, (t10, t11), (t20, t21))
