from __future__ import annotations

from typing import Any, Callable, TypeAlias, Iterable
from sympy import Abs, Expr, S, Derivative, Function as SymFunction, Min, Max, sympify, Add, Mul, Pow
from sympy.functions.elementary.miscellaneous import MinMaxBase
from sympy.physics.units import Dimension, Quantity as SymQuantity
from sympy.physics.units.systems.si import dimsys_SI

from .errors import UnitsError

ScalarValue: TypeAlias = Expr | float | int


class AnyDimension(Dimension):  # type: ignore[misc]
    # pylint: disable-next=signature-differs
    def __new__(cls) -> AnyDimension:
        return super().__new__(cls, "any_dimension")  # type: ignore[no-any-return]

    def _eval_nseries(self, _x: Any, _n: Any, _logx: Any, _cdir: Any) -> Any:
        pass


any_dimension = AnyDimension()


def _is_any_dimension(factor: Expr) -> bool:
    return factor in (S.Zero, S.Infinity, S.NegativeInfinity, S.NaN)


def _is_number(value: Any) -> bool:
    try:
        complex(value)
    except (TypeError, ValueError):
        return False

    return True


# pylint: disable-next=too-many-statements
def collect_quantity_factor_and_dimension(expr: Expr) -> tuple[Expr, Dimension]:
    """
    Returns tuple with scale factor expression and dimension expression. Designed to be
    used during the instantiation of the `symplyphysics.Quantity` class.

    Raises:
        ValueError: If the dimensions of the sub-expressions don't match.
        sympy.SympifyError: If ``expr`` cannot be converted to a value used by `sympy`.
    """

    def _collect_quantity(expr: SymQuantity) -> tuple[Expr, Dimension]:
        return (expr.scale_factor, expr.dimension)

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

        if _is_any_dimension(factor):
            return (factor, dimensionless)

        return (factor, dim * arg_dim)

    def _collect_pow(expr: Pow) -> tuple[Expr, Dimension]:
        (base_factor, base_dim) = collect_quantity_factor_and_dimension(expr.base)
        (exp_factor, exp_dim) = collect_quantity_factor_and_dimension(expr.exp)

        if _is_any_dimension(exp_factor) or dimsys_SI.is_dimensionless(exp_dim):
            return (base_factor**exp_factor, base_dim**exp_factor)

        raise ValueError(f"Dimension of '{expr.exp}' is {exp_dim}, but it should be dimensionless")

    @_elementwise_wrapper
    def _collect_add(factor: Expr, dim: Dimension, arg: Expr) -> tuple[Expr, Dimension]:
        arg_factor, arg_dim = collect_quantity_factor_and_dimension(arg)

        if _is_any_dimension(factor):
            dim = arg_dim
        elif _is_any_dimension(arg_factor):
            arg_dim = dim

        if not dimsys_SI.equivalent_dims(dim, arg_dim):
            raise ValueError(f"Dimension of '{arg}' is {arg_dim}, but it should be {dim}")

        return (factor + arg_factor, dim)

    def _collect_abs(expr: Abs) -> tuple[Expr, Dimension]:
        arg_factor, arg_dim = collect_quantity_factor_and_dimension(expr.args[0])
        return (Abs(arg_factor), arg_dim)

    def _collect_min_max(expr: MinMaxBase) -> tuple[Expr, Dimension]:
        cls = type(expr)

        @_elementwise_wrapper
        def collect(factor: Expr, dim: Dimension, arg: Expr) -> tuple[Expr, Dimension]:
            arg_factor, arg_dim = collect_quantity_factor_and_dimension(arg)

            if _is_any_dimension(factor):
                dim = arg_dim
            elif _is_any_dimension(arg_factor):
                arg_dim = dim

            if not dimsys_SI.equivalent_dims(dim, arg_dim):
                raise ValueError(f"Dimension of '{arg}' is {arg_dim}, but it should be {dim}")

            return (cls(factor, arg_factor), dim)

        return collect(expr)

    def _collect_function(expr: SymFunction) -> tuple[Expr, Dimension]:
        factors: list[Expr] = []

        for arg in expr.args:
            (arg_factor, arg_dim) = collect_quantity_factor_and_dimension(arg)

            # only functions with dimensionless arguments are supported
            if _is_any_dimension(arg_factor) or dimsys_SI.is_dimensionless(arg_dim):
                factors.append(arg_factor)
                continue

            raise ValueError(f"Dimension of '{arg}' is {arg_dim}, but it should be dimensionless")

        factor = expr.func(*(f for f in factors))
        return (factor, dimensionless)

    def _unsupported_derivative(expr: Derivative) -> tuple[Expr, Dimension]:
        raise ValueError(f"'{expr}' should not contain unevaluated Derivative")

    def _collect_default(expr: Expr) -> tuple[Expr, Dimension]:
        if not _is_number(expr):
            raise ValueError(f"'{expr}' should be an expression made of numbers or quantities.")

        return expr, dimensionless

    expr = sympify(expr)

    cases: dict[type, Callable[[Expr], tuple[Expr, Dimension]]] = {
        SymQuantity: _collect_quantity,
        Mul: _collect_mul,
        Pow: _collect_pow,
        Add: _collect_add,
        Abs: _collect_abs,
        MinMaxBase: _collect_min_max,
        Derivative: _unsupported_derivative,
        SymFunction: _collect_function,
    }

    for type_, collector in cases.items():
        if isinstance(expr, type_):
            return collector(expr)

    return _collect_default(expr)


def collect_expression_and_dimension(expr: Expr) -> tuple[Expr, Dimension]:
    """
    Returns the simplified representation and the dimension of the given expression. Unlike
    `collect_quantity_factor_and_dimension` it supports derivatives, symbols, and functions with
    dimensionful arguments.

    Raises:
        ValueError: if a sub-expression is dimensionful instead of dimensionless.
        UnitsError: if the dimensions of sub-expressions don't match.
    """

    from .symbols.quantities import Quantity

    def _split_numeric_and_symbolic(
        expr: Expr,) -> tuple[list[Expr], list[SymQuantity], list[tuple[Expr, Dimension]]]:
        nums: list[Expr] = []
        qtys: list[SymQuantity] = []
        syms: list[tuple[Expr, Dimension]] = []

        for arg in expr.args:
            if isinstance(arg, SymQuantity):
                qtys.append(arg)
            elif _is_number(arg):
                nums.append(arg)
            else:
                syms.append(collect_expression_and_dimension(arg))

        return nums, qtys, syms

    def _collect_mul(expr: Mul) -> tuple[Expr, Dimension]:
        nums, qtys, syms = _split_numeric_and_symbolic(expr)

        qty_factor = S.One
        for num in nums:
            qty_factor *= num

        qty_dim = dimensionless
        for qty in qtys:
            factor = qty.scale_factor
            qty_factor *= factor

            if not _is_any_dimension(factor):
                qty_dim *= qty.dimension

        if _is_any_dimension(qty_factor):
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

        if not _is_any_dimension(exp_expr) and not dimsys_SI.is_dimensionless(exp_dim):
            raise ValueError(
                f"Dimension of '{expr.exp}' is {exp_dim}, but it should be dimensionless")

        base_expr, base_dim = collect_expression_and_dimension(expr.base)

        expr_ = base_expr**exp_expr
        dim = base_dim**exp_expr

        return expr_, dim

    def _collect_unique_dimension(
        nums: Iterable[Expr],
        qtys: Iterable[SymQuantity],
        syms: Iterable[tuple[Expr, Dimension]],
    ) -> Dimension:
        dim = None

        if not all(_is_any_dimension(num) for num in nums):
            dim = dimensionless

        for qty in qtys:
            if dim is None:
                dim = qty.dimension
                continue

            if _is_any_dimension(qty.scale_factor):
                continue

            elif not dimsys_SI.equivalent_dims(dim, qty.dimension):
                raise UnitsError(f"The dimension of {qty} is {qty.dimension}, expected {dim}")

        for sym_expr, sym_dim in syms:
            if dim is None:
                dim = sym_dim
                continue

            if _is_any_dimension(sym_expr):
                continue

            elif not dimsys_SI.equivalent_dims(dim, sym_dim):
                raise UnitsError(f"The dimension of '{sym_expr}' is {sym_dim}, expected {dim}")

        # edge case when both `qtys` and `syms` are empty and all `nums` are of any dimension
        if dim is None:
            dim = dimensionless

        return dim

    def _collect_add(expr: Add) -> tuple[Expr, Dimension]:
        nums, qtys, syms = _split_numeric_and_symbolic(expr)
        dim = _collect_unique_dimension(nums, qtys, syms)

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

        nums, qtys, syms = _split_numeric_and_symbolic(expr)
        dim = _collect_unique_dimension(nums, qtys, syms)

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

    expr = sympify(expr)

    # early return that works for `sympy.Quantity`, `SymbolNew`, `SymbolIndexedNew`, and `FunctionNew`
    if hasattr(expr, "dimension"):
        return expr, getattr(expr, "dimension")

    cases: dict[type, Callable[[Any], tuple[Expr, Dimension]]] = {
        Mul: _collect_mul,
        Pow: _collect_pow,
        Add: _collect_add,
        Abs: _collect_abs,
        Min: _collect_min_max,
        Max: _collect_min_max,
        Derivative: _collect_derivative,
        SymFunction: _collect_function,
    }

    for type_, collector in cases.items():
        if isinstance(expr, type_):
            return collector(expr)

    return (expr, dimensionless)


def assert_equivalent_dimension(
    arg: SymQuantity | ScalarValue | Dimension,
    param_name: str,
    func_name: str,
    expected_unit: SymQuantity | Dimension,
) -> None:
    """
    Asserts if the dimension of the argument matches the provided unit.

    Args:
        arg: Number, quantity, expression made of numbers and quantities, or dimension.
        param_name: Name of the parameter of the calling function.
        func_name: Name of the calling function.
        expected_unit: Expression or dimension which `arg` is compared to.

    Raises:
        TypeError: If `arg` is a number, but `expected_unit` is not dimensionless.
        UnitsError: If the dimensions don't match otherwise, or when the scale factor of `arg` is not a number.
    """

    if not isinstance(expected_unit, Dimension):
        expected_scale_factor, expected_unit = collect_quantity_factor_and_dimension(expected_unit)

        if _is_any_dimension(expected_scale_factor) or isinstance(expected_unit, AnyDimension):
            return

    # HACK: this allows to treat angle type as dimensionless
    expected_unit = expected_unit.subs("angle", S.One)

    if not isinstance(arg, Dimension):
        (scale_factor, arg) = collect_quantity_factor_and_dimension(arg)

        if not _is_number(scale_factor):
            # NOTE: this should probably raise `ValueError` or `TypeError`
            raise UnitsError(f"Argument '{param_name}' to function '{func_name}' should "
                f"not contain free symbols: '{scale_factor}'")

        if _is_any_dimension(scale_factor) or isinstance(arg, AnyDimension):
            return

    # HACK: this allows to treat angle type as dimensionless
    arg = arg.subs("angle", S.One)

    if dimsys_SI.is_dimensionless(arg) and not dimsys_SI.is_dimensionless(expected_unit):
        # NOTE: this should probably be `UnitsError`
        raise TypeError(f"Argument '{param_name}' to function '{func_name}'"
            f" is Number but '{expected_unit}' is not dimensionless")

    if not dimsys_SI.equivalent_dims(arg, expected_unit):
        raise UnitsError(f"Argument '{param_name}' to function '{func_name}' must "
            f"be in units equivalent to '{expected_unit.name}', got {arg.name}")


dimensionless = Dimension(S.One)


def print_dimension(dimension: Dimension) -> str:
    return "dimensionless" if dimsys_SI.is_dimensionless(dimension) else str(dimension.name)


__all__ = [
    # re-exports
    "Dimension",

    # locals
    "ScalarValue",
    "AnyDimension",
    "any_dimension",
    "collect_quantity_factor_and_dimension",
    "collect_expression_and_dimension",
    "assert_equivalent_dimension",
    "dimensionless",
    "print_dimension",
]
