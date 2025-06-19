from __future__ import annotations

from typing import Optional, Any

from sympy import Expr, Symbol as SymSymbol, Integral as SymIntegral, diff
from sympy.physics import units

from symplyphysics.core.symbols.symbols import Symbol

from ..miscellaneous import evaluate_or_global_fallback

from ..vectors import VectorSymbol, VectorNorm
from ..coordinate_systems.curve import Curve
from ..coordinate_systems.vector import CoordinateVector

INFINITESIMAL_DISPLACEMENT = VectorSymbol("dr", units.length, display_latex="d \\! \\vec r")
"""To be used when the integrand is a vector field."""

INFINITESIMAL_ARC_LENGTH = Symbol("ds", units.length, display_latex="d \\! s", positive=True)
"""To be used when the integrand is a scalar field."""


class LineIntegral(Expr):  # pylint: disable=too-few-public-methods
    """
    A **line integral**, or **path integral**, is an integral where the integrand is evaluated along a
    curve.

    Line integral `formula <https://en.wikipedia.org/wiki/Line_integral#Line_integral_of_a_scalar_field>`__
    when the integrand is a scalar field.

    Line integral `formula <https://en.wikipedia.org/wiki/Line_integral#Line_integral_of_a_vector_field>`__
    when the integrand is a vector field.

    More formulas can be found `here <https://en.wikipedia.org/wiki/Curvilinear_coordinates#Integration_2>`__,
    see the first row of the table.

    ..
        NOTE: we do not distinguish yet between line integrals over open and closed curves
    """

    def __new__(
        cls,
        expr: Expr,
        curve: Curve | SymSymbol,
        bounds: Optional[tuple[Any, Any]] = None,
        *,
        evaluate: Optional[bool] = None,
    ) -> Expr:
        evaluate = evaluate_or_global_fallback(evaluate)

        if not evaluate or isinstance(curve, SymSymbol):
            if not bounds:
                return super().__new__(cls, expr, curve)  # pylint: disable=too-many-function-args

            return super().__new__(cls, expr, curve, bounds)  # pylint: disable=too-many-function-args

        if bounds is None:
            raise ValueError("Define the curve parameter bounds for the line integral")

        t = curve.parameter

        system = curve.parametrization.system
        q1, q2, q3 = curve.parametrization.coordinates.values()

        expr = expr.subs(system.base_scalar_subs((q1, q2, q3)))

        h1, h2, h3 = system.lame_coefficients((q1, q2, q3))

        dr = CoordinateVector(
            [h1 * diff(q1, t), h2 * diff(q2, t), h3 * diff(q3, t)],
            system,
            curve.parametrization,
        )

        ds = VectorNorm(dr)

        expr = expr.subs({
            INFINITESIMAL_DISPLACEMENT: dr,
            INFINITESIMAL_ARC_LENGTH: ds,
        })

        t0, t1 = bounds
        return SymIntegral(expr, (t, t0, t1))


__all__ = [
    "INFINITESIMAL_DISPLACEMENT",
    "INFINITESIMAL_ARC_LENGTH",
    "LineIntegral",
]
