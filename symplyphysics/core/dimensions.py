from __future__ import annotations

from typing import Any, Callable, TypeAlias, Iterable
from sympy import Abs, Expr, S, Derivative, Function as SymFunction, Min, Max, sympify, Add, Mul, Pow
from sympy.functions.elementary.miscellaneous import MinMaxBase
from sympy.physics.units import Dimension, Quantity as SymQuantity
from sympy.physics.units.systems.si import dimsys_SI

from .errors import UnitsError

ScalarValue: TypeAlias = Expr | float


class AnyDimension(Dimension):  # type: ignore[misc]
    # pylint: disable-next=signature-differs
    def __new__(cls) -> AnyDimension:
        return super().__new__(cls, "any_dimension")  # type: ignore[no-any-return]

    def _eval_nseries(self, _x: Any, _n: Any, _logx: Any, _cdir: Any) -> Any:
        pass


any_dimension = AnyDimension()


def _is_number(value: Any) -> bool:
    try:
        complex(value)
    except (TypeError, ValueError):
        return False

    return True


# pylint: disable-next=too-many-statements
def collect_quantity_factor_and_dimension(expr: Expr) -> tuple[Expr, Dimension]:
    """
    Return tuple with scale factor expression and dimension expression. Designed to be
    used during the instantiation of the `symplyphysics.Quantity` class.
    """

    def _is_any_dimension(factor: Expr) -> bool:
        return factor in (S.Zero, S.Infinity, S.NegativeInfinity)

    def _collect_quantity(expr: SymQuantity) -> tuple[Expr, Dimension]:
        return (expr.scale_factor, expr.dimension)

    def _collect_binary_wrapper(
        inner: Callable[[Expr, Dimension, Expr], tuple[Expr, Dimension]]
    ) -> Callable[[Expr], tuple[Expr, Dimension]]:

        def outer(expr: Expr) -> tuple[Expr, Dimension]:
            factor, dim = collect_quantity_factor_and_dimension(expr.args[0])
            for arg in expr.args[1:]:
                factor, dim = inner(factor, dim, arg)
            return (factor, dim)

        return outer

    @_collect_binary_wrapper
    def _collect_mul(factor: Expr, dim: Dimension, arg: Expr) -> tuple[Expr, Dimension]:
        arg_factor, arg_dim = collect_quantity_factor_and_dimension(arg)

        return (factor * arg_factor, dim * arg_dim)

    def _collect_pow(expr: Pow) -> tuple[Expr, Dimension]:
        (base_factor, base_dim) = collect_quantity_factor_and_dimension(expr.base)
        (exp_factor, exp_dim) = collect_quantity_factor_and_dimension(expr.exp)

        if _is_any_dimension(exp_factor) or dimsys_SI.is_dimensionless(exp_dim):
            return (base_factor**exp_factor, base_dim**exp_factor)

        raise ValueError(f"Dimension of '{expr.exp}' is {exp_dim}, but it should be dimensionless")

    @_collect_binary_wrapper
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

        @_collect_binary_wrapper
        def wrapped(factor: Expr, dim: Dimension, arg: Expr) -> tuple[Expr, Dimension]:
            arg_factor, arg_dim = collect_quantity_factor_and_dimension(arg)

            if _is_any_dimension(factor):
                dim = arg_dim
            elif _is_any_dimension(arg_factor):
                arg_dim = dim

            if not dimsys_SI.equivalent_dims(dim, arg_dim):
                raise ValueError(f"Dimension of '{arg}' is {arg_dim}, but it should be {dim}")

            return (cls(factor, arg_factor), dim)

        return wrapped(expr)

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


def collect_dimension(expr: Expr) -> Dimension:
    """
    Returns the dimension of the given expression. Unlike `collect_quantity_factor_and_dimension`
    it supports derivatives, symbols, and functions with dimensionful arguments.
    """

    def _collect_mul(expr: Mul) -> Dimension:
        dim = dimensionless
        for arg in expr.args:
            if arg == S.Zero:
                return dimensionless

            dim *= collect_dimension(arg)
        return dim

    def _collect_pow(expr: Pow) -> Dimension:
        exp_dim = collect_dimension(expr.exp)
        if not dimsys_SI.is_dimensionless(exp_dim):
            raise ValueError(
                f"Dimension of '{expr.exp}' is {exp_dim}, but it should be dimensionless")

        if expr.base == S.Zero or expr.exp == S.Zero:  # pylint: disable=consider-using-in
            return Dimension(S.One)

        base_dim = collect_dimension(expr.base)

        # NOTE: this works fine when `expr.exp` is just a number or a `sympy.Symbol`
        #       but (see line 206) `dimsys_SI.is_dimensionless` breaks if it's a complex
        #       expression (e.g. `time * temporal_frequency`), although the dimension
        #       itself evaluates as usual
        return base_dim**expr.exp

    def _collect_unique_dimension(args: Iterable[Expr]) -> Dimension:
        dims = set(collect_dimension(arg) for arg in args)
        if len(dims) > 1:
            raise ValueError(f"Arguments must have the same dimension, got {tuple(dims)}")
        (dim,) = dims
        return dim

    def _collect_add(expr: Add) -> Dimension:
        return _collect_unique_dimension(expr.args)

    def _collect_abs(expr: Abs) -> Dimension:
        return collect_dimension(expr.args[0])

    def _collect_min(expr: Min) -> Dimension:
        return _collect_unique_dimension(expr.args)

    def _collect_derivative(expr: Derivative) -> Dimension:
        func, *args = expr.args
        dim = collect_dimension(func.func)
        for arg, n in args:
            arg_dim = collect_dimension(arg)
            dim /= arg_dim**n
        return dim

    # early return that works for `sympy.Quantity`, `SymbolNew`, `SymbolIndexedNew`, and `FunctionNew`
    if hasattr(expr, "dimension"):
        return getattr(expr, "dimension")

    cases: dict[type, Callable[[Any], tuple[Expr, Dimension]]] = {
        Mul: _collect_mul,
        Pow: _collect_pow,
        Add: _collect_add,
        Abs: _collect_abs,
        Min: _collect_min,
        Derivative: _collect_derivative,
    }

    for k, v in cases.items():
        if isinstance(expr, k):
            dim = v(expr)

            if dimsys_SI.is_dimensionless(dim):
                return dimensionless

            return dim

    return dimensionless


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

    def _is_any_dimension(scale_factor: Expr, dimension: Dimension) -> bool:
        """`Zero` and `Infinity` can be of any dimension, as well as `AnyDimension` instances."""

        return (scale_factor in (S.Zero, S.Infinity, S.NegativeInfinity) or
            isinstance(dimension, AnyDimension))

    if not isinstance(expected_unit, Dimension):
        expected_scale_factor, expected_unit = collect_quantity_factor_and_dimension(expected_unit)

        if _is_any_dimension(expected_scale_factor, expected_unit):
            return

    # HACK: this allows to treat angle type as dimensionless
    expected_unit = expected_unit.subs("angle", S.One)

    if not isinstance(arg, Dimension):
        (scale_factor, arg) = collect_quantity_factor_and_dimension(sympify(arg))

        if not _is_number(scale_factor):
            # NOTE: this should probably raise `ValueError` or `TypeError`
            raise UnitsError(f"Argument '{param_name}' to function '{func_name}' should "
                f"not contain free symbols: '{scale_factor}'")

        if _is_any_dimension(scale_factor, arg):
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
    "collect_dimension",
    "assert_equivalent_dimension",
    "dimensionless",
    "print_dimension",
]
