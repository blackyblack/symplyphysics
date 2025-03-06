from typing import Callable, SupportsFloat
from sympy import Expr, Pow, Derivative, Abs, Mul, Add, Function as SymFunction, sympify
from sympy.functions.elementary.miscellaneous import MinMaxBase
from sympy.physics.units import Quantity as SymQuantity, Dimension
from sympy.physics.units.prefixes import Prefix
from sympy.physics.units.systems.si import dimsys_SI

from .miscellaneous import is_any_dimension, is_number, dimensionless


def _collect_quantity(expr: SymQuantity) -> tuple[Expr, Dimension]:
    return (expr.scale_factor, expr.dimension)


def _collect_prefix(expr: Prefix) -> tuple[Expr, Dimension]:
    return (expr.scale_factor, dimensionless)


def _elementwise_wrapper(
    inner: Callable[[Expr, Dimension, Expr], tuple[Expr, Dimension]]
) -> Callable[[Expr], tuple[Expr, Dimension]]:

    def outer(expr: Expr) -> tuple[Expr, Dimension]:
        factor, dim = collect_quantity_factor_and_dimension(expr.args[0])
        for arg in expr.args[1:]:
            factor, dim = inner(factor, dim, arg)
        return (factor, dim)

    return outer


@_elementwise_wrapper
def _collect_mul(factor: Expr, dim: Dimension, arg: Expr) -> tuple[Expr, Dimension]:
    arg_factor, arg_dim = collect_quantity_factor_and_dimension(arg)

    factor *= arg_factor

    if is_any_dimension(factor):
        return (factor, dimensionless)

    return (factor, dim * arg_dim)


def _collect_pow(expr: Pow) -> tuple[Expr, Dimension]:
    (base_factor, base_dim) = collect_quantity_factor_and_dimension(expr.base)
    (exp_factor, exp_dim) = collect_quantity_factor_and_dimension(expr.exp)

    if is_any_dimension(exp_factor) or dimsys_SI.is_dimensionless(exp_dim):
        return (base_factor**exp_factor, base_dim**exp_factor)

    raise ValueError(f"Dimension of '{expr.exp}' is {exp_dim}, but it should be dimensionless")


@_elementwise_wrapper
def _collect_add(factor: Expr, dim: Dimension, arg: Expr) -> tuple[Expr, Dimension]:
    arg_factor, arg_dim = collect_quantity_factor_and_dimension(arg)

    if is_any_dimension(factor):
        dim = arg_dim
    elif is_any_dimension(arg_factor):
        arg_dim = dim

    if not dimsys_SI.equivalent_dims(dim, arg_dim):
        raise ValueError(f"Dimension of '{arg}' is {arg_dim}, but it should be {dim}")

    return (factor + arg_factor, dim)


def _collect_abs(expr: Abs) -> tuple[Expr, Dimension]:
    arg_factor, arg_dim = collect_quantity_factor_and_dimension(expr.args[0])
    return (Abs(arg_factor), arg_dim)


def _collect_min_max(expr: MinMaxBase) -> tuple[Expr, Dimension]:
    cls = type(expr)

    def collect(factor: Expr, dim: Dimension, arg: Expr) -> tuple[Expr, Dimension]:
        arg_factor, arg_dim = collect_quantity_factor_and_dimension(arg)

        if is_any_dimension(factor):
            dim = arg_dim
        elif is_any_dimension(arg_factor):
            arg_dim = dim

        if not dimsys_SI.equivalent_dims(dim, arg_dim):
            raise ValueError(f"Dimension of '{arg}' is {arg_dim}, but it should be {dim}")

        return (cls(factor, arg_factor), dim)

    return _elementwise_wrapper(collect)(expr)


def _collect_function(expr: SymFunction) -> tuple[Expr, Dimension]:
    factors: list[Expr] = []

    for arg in expr.args:
        (arg_factor, arg_dim) = collect_quantity_factor_and_dimension(arg)

        # only functions with dimensionless arguments are supported
        if is_any_dimension(arg_factor) or dimsys_SI.is_dimensionless(arg_dim):
            factors.append(arg_factor)
            continue

        raise ValueError(f"Dimension of '{arg}' is {arg_dim}, but it should be dimensionless")

    factor = expr.func(*(f for f in factors))
    return (factor, dimensionless)


def _unsupported_derivative(expr: Derivative) -> tuple[Expr, Dimension]:
    raise ValueError(f"'{expr}' should not contain unevaluated Derivative")


def _collect_default(expr: Expr) -> tuple[Expr, Dimension]:
    if not is_number(expr):
        raise ValueError(f"'{expr}' should be an expression made of numbers or quantities.")

    return expr, dimensionless


_cases: dict[type, Callable[[Expr], tuple[Expr, Dimension]]] = {
    SymQuantity: _collect_quantity,
    Prefix: _collect_prefix,
    Mul: _collect_mul,
    Pow: _collect_pow,
    Add: _collect_add,
    Abs: _collect_abs,
    MinMaxBase: _collect_min_max,
    Derivative: _unsupported_derivative,
    SymFunction: _collect_function,
}


def collect_quantity_factor_and_dimension(expr: SupportsFloat) -> tuple[Expr, Dimension]:
    """
    Returns tuple with scale factor expression and dimension expression. Designed to be
    used during the instantiation of the `symplyphysics.Quantity` class.

    Raises:
        ValueError: If the dimensions of the sub-expressions don't match.
        sympy.SympifyError: If ``expr`` cannot be converted to a value used by `sympy`.
    """

    expr = sympify(expr)  # no `strict` because we want to allow sympy functions here too

    for type_, collector in _cases.items():
        if isinstance(expr, type_):
            return collector(expr)

    return _collect_default(expr)


__all__ = ["collect_quantity_factor_and_dimension"]
