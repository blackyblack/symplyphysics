from __future__ import annotations

from typing import Iterable, Callable, Any, SupportsFloat
from sympy import Expr, S, Mul, Add, Pow, Abs, Min, Max, Derivative, Function as SymFunction, sympify, Basic
from sympy.functions.elementary.miscellaneous import MinMaxBase
from sympy.physics import units
from sympy.physics.units import Dimension, Quantity as SymQuantity
from sympy.physics.units.systems.si import dimsys_SI

from ..experimental.vectors import (VectorDot, VectorCross, VectorNorm, VectorMixedProduct,
    AppliedVectorFunction)
from ..experimental.coordinate_systems import CoordinateScalar, CoordinateVector
from ..experimental.operators import VectorGradient, VectorDivergence, VectorCurl, VectorLaplacian

from ..errors import UnitsError
from ..symbols.quantities import Quantity
from .miscellaneous import is_number, is_any_dimension, dimensionless


def split_numeric_and_symbolic(
    args: Iterable[Expr],) -> tuple[list[Expr], list[SymQuantity], list[tuple[Expr, Dimension]]]:
    nums: list[Expr] = []
    qtys: list[SymQuantity] = []
    syms: list[tuple[Expr, Dimension]] = []

    for arg in args:
        if isinstance(arg, SymQuantity):
            qtys.append(arg)
        elif is_number(arg):
            nums.append(arg)
        else:
            syms.append(collect_expression_and_dimension(arg))

    return nums, qtys, syms


def _collect_mul(expr: Mul) -> tuple[Expr, Dimension]:
    nums, qtys, syms = split_numeric_and_symbolic(expr.args)

    qty_factor = S.One
    for num in nums:
        qty_factor *= num

    qty_dim = dimensionless
    for qty in qtys:
        factor = qty.scale_factor
        qty_factor *= factor

        if not is_any_dimension(factor):
            qty_dim *= qty.dimension

    if is_any_dimension(qty_factor):
        return qty_factor, dimensionless

    if dimsys_SI.is_dimensionless(qty_dim):
        expr_ = qty_factor
    else:
        expr_ = Quantity(qty_factor, dimension=qty_dim)

    dim = qty_dim

    for sym_expr, sym_dim in syms:
        expr_ *= sym_expr
        dim *= sym_dim

    return expr_, dim


def _collect_pow(expr: Pow) -> tuple[Expr, Dimension]:
    exp_expr, exp_dim = collect_expression_and_dimension(expr.exp)

    if not is_any_dimension(exp_expr) and not dimsys_SI.is_dimensionless(exp_dim):
        raise ValueError(f"Dimension of '{expr.exp}' is {exp_dim}, but it should be dimensionless")

    base_expr, base_dim = collect_expression_and_dimension(expr.base)

    expr_ = base_expr**exp_expr
    dim = base_dim**exp_expr

    return expr_, dim


def collect_unique_dimension(
    nums: Iterable[Expr],
    qtys: Iterable[SymQuantity],
    syms: Iterable[tuple[Expr, Dimension]],
) -> Dimension:
    dim = None

    if not all(is_any_dimension(num) for num in nums):
        dim = dimensionless

    for qty in qtys:
        if dim is None:
            dim = qty.dimension
            continue

        if is_any_dimension(qty.scale_factor):
            continue

        if not dimsys_SI.equivalent_dims(dim, qty.dimension):
            raise UnitsError(f"The dimension of {qty} is {qty.dimension}, expected {dim}")

    for sym_expr, sym_dim in syms:
        if dim is None:
            dim = sym_dim
            continue

        if is_any_dimension(sym_expr):
            continue

        if not dimsys_SI.equivalent_dims(dim, sym_dim):
            raise UnitsError(f"The dimension of '{sym_expr}' is {sym_dim}, expected {dim}")

    # edge case when both `qtys` and `syms` are empty and all `nums` are of any dimension
    if dim is None:
        dim = dimensionless

    return dim


def _collect_add(expr: Add) -> tuple[Expr, Dimension]:
    nums, qtys, syms = split_numeric_and_symbolic(expr.args)
    dim = collect_unique_dimension(nums, qtys, syms)

    qty_sum = sum(nums, start=S.Zero) + sum((qty.scale_factor for qty in qtys), start=S.Zero)

    if not dimsys_SI.is_dimensionless(dim):
        qty_sum = Quantity(qty_sum, dimension=dim)

    sym_sum = sum((sym_expr for (sym_expr, _) in syms), start=S.Zero)

    return qty_sum + sym_sum, dim


def _collect_abs(expr: Abs) -> tuple[Expr, Dimension]:
    expr_, dim = collect_expression_and_dimension(expr.args[0])
    return Abs(expr_), dim


def _collect_min_max(expr: MinMaxBase) -> tuple[Expr, Dimension]:
    cls = type(expr)

    nums, qtys, syms = split_numeric_and_symbolic(expr.args)
    dim = collect_unique_dimension(nums, qtys, syms)

    if dimsys_SI.is_dimensionless(dim):
        expr_ = cls(
            *nums,
            *(qty.scale_factor for qty in qtys),
            *(sym_expr for (sym_expr, _) in syms),
        )
        return expr_, dim

    minmax_num = cls(*nums)
    minmax_qty_factor = cls(minmax_num, *(qty.scale_factor for qty in qtys))
    minmax_qty = Quantity(minmax_qty_factor, dimension=dim)

    expr_ = cls(minmax_qty, *(sym_expr for (sym_expr, _) in syms))

    return expr_, dim


def _collect_function(expr: SymFunction) -> tuple[Expr, Dimension]:
    func = expr.func
    dim = getattr(func, "dimension", dimensionless)

    factors = []
    for arg in expr.args:
        arg_factor, _ = collect_expression_and_dimension(arg)
        factors.append(arg_factor)

    expr_ = func(*factors)
    return expr_, dim


def _collect_derivative(expr: Derivative) -> tuple[Expr, Dimension]:
    func, *args = expr.args
    _, dim = collect_expression_and_dimension(func.func)

    expr_ = func
    for arg, n in args:
        arg_expr, arg_dim = collect_expression_and_dimension(arg)
        dim /= arg_dim**n
        expr_ = expr_.diff((arg_expr, n))
    return expr_, dim


def _collect_vector_dot(expr: VectorDot) -> tuple[Expr, Dimension]:
    if 0 in expr.args:
        return S.Zero, dimensionless

    lhs, rhs = expr.args

    _, lhs_dim = collect_expression_and_dimension(lhs)
    _, rhs_dim = collect_expression_and_dimension(rhs)

    return expr, lhs_dim * rhs_dim


def _collect_vector_cross(expr: VectorCross) -> tuple[Expr, Dimension]:
    if 0 in expr.args:
        return S.Zero, dimensionless

    lhs, rhs = expr.args

    _, lhs_dim = collect_expression_and_dimension(lhs)
    _, rhs_dim = collect_expression_and_dimension(rhs)

    return expr, lhs_dim * rhs_dim


def _collect_vector_mixed_product(expr: VectorMixedProduct) -> tuple[Expr, Dimension]:
    if 0 in expr.args:
        return S.Zero, dimensionless

    a, b, c = expr.args

    _, a_dim = collect_expression_and_dimension(a)
    _, b_dim = collect_expression_and_dimension(b)
    _, c_dim = collect_expression_and_dimension(c)

    return expr, a_dim * b_dim * c_dim


def _collect_vector_norm(expr: VectorNorm) -> tuple[Expr, Dimension]:
    if 0 in expr.args:
        return S.Zero, dimensionless

    _, dim = collect_expression_and_dimension(expr.args[0])

    return expr, dim


def _collect_applied_vector_function(expr: AppliedVectorFunction) -> tuple[Expr, Dimension]:
    _, dim = collect_expression_and_dimension(expr.func)

    return expr, dim


def _collect_coordinate_scalar(expr: CoordinateScalar) -> tuple[Expr, Dimension]:
    _, dim = collect_expression_and_dimension(expr.scalar)

    return expr, dim


def _collect_coordinate_vector(expr: CoordinateVector) -> tuple[Expr, Dimension]:
    nums, qtys, syms = split_numeric_and_symbolic(expr.components)
    dim = collect_unique_dimension(nums, qtys, syms)

    return expr, dim


def _collect_grad_div_curl(
        expr: VectorGradient | VectorDivergence | VectorCurl) -> tuple[Expr, Dimension]:
    _, dim = collect_expression_and_dimension(expr.args[0])

    return expr, dim / units.length


def _collect_laplacian(expr: VectorLaplacian) -> tuple[Expr, Dimension]:
    _, dim = collect_expression_and_dimension(expr.args[0])

    return expr, dim / units.length**2


def _collect_dimension(expr: Dimension) -> tuple[Expr, Dimension]:
    raise TypeError(f"{expr} cannot appear in an expression context")


_cases: dict[type, Callable[[Any], tuple[Expr, Dimension]]] = {
    Mul: _collect_mul,
    Pow: _collect_pow,
    Add: _collect_add,
    Abs: _collect_abs,
    Min: _collect_min_max,
    Max: _collect_min_max,
    Derivative: _collect_derivative,
    SymFunction: _collect_function,
    VectorDot: _collect_vector_dot,
    VectorCross: _collect_vector_cross,
    VectorMixedProduct: _collect_vector_mixed_product,
    VectorNorm: _collect_vector_norm,
    AppliedVectorFunction: _collect_applied_vector_function,
    CoordinateScalar: _collect_coordinate_scalar,
    CoordinateVector: _collect_coordinate_vector,
    VectorGradient: _collect_grad_div_curl,
    VectorDivergence: _collect_grad_div_curl,
    VectorCurl: _collect_grad_div_curl,
    VectorLaplacian: _collect_laplacian,
    Dimension: _collect_dimension,
}


def collect_expression_and_dimension(expr: SupportsFloat) -> tuple[Expr, Dimension]:
    """
    Returns the simplified representation and the dimension of the given expression. Unlike
    `collect_quantity_factor_and_dimension` it supports derivatives, symbols, and functions with
    dimensionful arguments.

    Raises:
        ValueError: if a sub-expression is dimensionful instead of dimensionless.
        UnitsError: if the dimensions of sub-expressions don't match.
    """

    expr = sympify(expr)  # no `strict` because we want to allow sympy functions here too

    # early return that works for `sympy.Quantity`, `SymbolNew`, `SymbolIndexedNew`, and `FunctionNew`
    if hasattr(expr, "dimension"):
        return expr, getattr(expr, "dimension")

    if isinstance(expr, Basic):
        expr = expr.doit()

    for type_, collector in _cases.items():
        if isinstance(expr, type_):
            return collector(expr)

    return (expr, dimensionless)


__all__ = [
    "collect_expression_and_dimension",
]
