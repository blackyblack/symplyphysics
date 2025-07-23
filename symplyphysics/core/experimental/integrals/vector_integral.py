from typing import Sequence

from sympy import Expr, Symbol as SymSymbol

from ..vectors import into_terms, is_vector_expr
from ..coordinate_systems import (combine_coordinate_vectors, CoordinateVector,
    CartesianCoordinateSystem, QuantityCoordinateVector)


def _check_symbols(symbols: Sequence[SymSymbol | tuple[SymSymbol, Expr, Expr]]) -> None:
    for symbol_or_tuple in symbols:
        if isinstance(symbol_or_tuple, SymSymbol):
            symbol = symbol_or_tuple
        else:
            symbol, _, _ = symbol_or_tuple

        if is_vector_expr(symbol):
            raise ValueError(f"Expected a scalar symbol, got vector {symbol}")


def integrate_coordinate_vectors(
    expr: Expr,
    *symbols: SymSymbol | tuple[SymSymbol, Expr, Expr],
) -> Expr:
    """
    Integrates the vector `expr` with respect to given `symbols`. The given `expr` must be composed
    of `CoordinateVector` instances defined in a `CartesianCoordinateSystem`. Note that a free
    vector term is not added to the resulting computation.
    """

    _check_symbols(symbols)

    result_sum = 0

    for term in into_terms(combine_coordinate_vectors(expr)):
        if not isinstance(term, CoordinateVector):
            raise TypeError(
                "Only `CoordinateVector` instances are supported for vector integration")

        # NOTE: technically the integration is possible if the point of vector application doesn't
        # change, but it's hard to represent this in code, so we disallow it altogether.
        if not isinstance(term.system, CartesianCoordinateSystem):
            raise TypeError("Expected a Cartesian vector; consider converting non-Cartesian "
                f"({type(term.system).__name__}) vectors to a Cartesian representation first")

        result_matrix = term.components.integrate(*symbols)

        if isinstance(term, QuantityCoordinateVector):
            try:
                result_term = QuantityCoordinateVector(result_matrix, term.system, term.point)
            except ValueError:
                result_term = CoordinateVector(result_matrix, term.system, term.point)
        else:
            result_term = CoordinateVector(result_matrix, term.system, term.point)

        result_sum += result_term

    return combine_coordinate_vectors(result_sum)
