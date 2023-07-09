from __future__ import annotations
from functools import partial
from typing import Any, Callable, Optional
from sympy import Expr, S, Derivative, Function as SymFunction, exp, Basic, sympify
from sympy.core.add import Add
from sympy.core.mul import Mul
from sympy.core.power import Pow
from sympy.physics.units import Dimension, Quantity as SymQuantity
from sympy.physics.units.systems.si import SI
from .symbols import DimensionSymbol, next_name


class Quantity(DimensionSymbol, SymQuantity):  # pylint: disable=too-many-ancestors

    def __new__(cls,
        _name: Basic,
        _abbrev: Optional[str] = None,
        _latex_repr: Optional[str] = None,
        _pretty_unicode_repr: Optional[str] = None,
        _pretty_ascii_repr: Optional[str] = None,
        _mathml_presentation_repr: Optional[str] = None,
        _is_prefixed: bool = False,
        **assumptions: Any) -> Quantity:
        name = next_name("QTY")
        obj = SymQuantity.__new__(cls, name, None, None, None, None, None, False, **assumptions)
        return obj

    def __init__(self, expr: Basic | float = S.One, *, dimension: Optional[Dimension] = None):
        (scale, dimension_) = collect_factor_and_dimension(sympify(expr))
        dimension = dimension_ if dimension is None else dimension
        super().__init__(self.name, dimension)
        SI.set_quantity_dimension(self, dimension)
        SI.set_quantity_scale_factor(self, scale)

    # This is required for integration to work properly
    @property
    def func(self) -> partial:
        return partial(Quantity.identity, self)

    def identity(self, *_args: Any) -> Quantity:
        return self


# HACK: copy of SI._collect_factor_and_dimension with fixed exp() dimension evaluation
def collect_factor_and_dimension(expr: Basic) -> tuple[Basic, Dimension]:
    """
    Return tuple with scale factor expression and dimension expression.
    """

    def _collect_quantity(expr: SymQuantity) -> tuple[Basic, Dimension]:
        return (expr.scale_factor, expr.dimension)

    def _collect_mul(expr: Mul) -> tuple[Basic, Dimension]:
        factor = S.One
        dimension = Dimension(S.One)
        for arg in expr.args:
            (arg_factor, arg_dim) = collect_factor_and_dimension(arg)
            factor *= arg_factor
            dimension *= arg_dim
        return (factor, dimension)

    def _collect_pow(expr: Pow) -> tuple[Basic, Dimension]:
        pow_expr: Expr = S.One
        (factor, dim) = collect_factor_and_dimension(expr.base)
        pow_expr *= factor
        (exp_factor, exp_dim) = collect_factor_and_dimension(expr.exp)
        if not SI.get_dimension_system().is_dimensionless(exp_dim):
            raise ValueError(
                f"Dimension of '{expr.exp}' is {exp_dim}, but it should be dimensionless")
        exp_dim = S.One
        return (pow_expr**exp_factor, dim**(exp_factor * exp_dim))

    def _collect_add(expr: Add) -> tuple[Basic, Dimension]:
        sum_expr: Expr = S.Zero
        (factor, dim) = collect_factor_and_dimension(expr.args[0])
        sum_expr += factor
        for addend in expr.args[1:]:
            (addend_factor, addend_dim) = collect_factor_and_dimension(addend)
            # automatically convert zero to the dimension of it's additives
            if dim != addend_dim and not SI.get_dimension_system().equivalent_dims(dim, addend_dim):
                if factor == S.Zero:
                    dim = addend_dim
                elif addend_factor == S.Zero:
                    addend_dim = dim
            if dim != addend_dim and not SI.get_dimension_system().equivalent_dims(dim, addend_dim):
                raise ValueError(f"Dimension of '{addend}' is {addend_dim}, but it should be {dim}")
            sum_expr += addend_factor
        return (sum_expr, dim)

    def _unsupported_derivative(expr: Derivative):
        raise ValueError(f"Dimension '{expr}' should not contain unevaluated Derivative")

    def _collect_function(expr: SymFunction) -> tuple[Basic, Dimension]:
        factors: list[Basic] = []
        for arg in expr.args:
            (f, d) = collect_factor_and_dimension(arg)
            # only functions with dimensionless arguments are supported
            if not SI.get_dimension_system().is_dimensionless(d):
                raise ValueError(f"Dimension of '{arg}' is {d}, but it should be dimensionless")
            factors.append(f)
        ret = expr.func(*(f for f in factors))
        return (ret, dimensionless)

    def _collect_dimension(expr: Dimension) -> tuple[Basic, Dimension]:
        return (S.One, expr)

    cases: dict[type, Callable[[Any], tuple[Basic, Dimension]]] = {
        SymQuantity: _collect_quantity,
        Mul: _collect_mul,
        Pow: _collect_pow,
        Add: _collect_add,
        Derivative: _unsupported_derivative,
        SymFunction: _collect_function,
        Dimension: _collect_dimension,
    }

    for k, v in cases.items():
        if isinstance(expr, k):
            return v(expr)
    return (expr, dimensionless)


def expr_to_quantity(expr: Expr) -> Quantity:
    quantity_scale = collect_factor_and_dimension(expr)
    dimension = quantity_scale[1]
    # HACK: this allows to treat angle type as dimensionless
    dimension = dimension.subs("angle", S.One)
    if isinstance(dimension, exp):
        dimension_args = dimension.nargs.args[0]
        # power of the exponent should be dimensionless
        assert dimension_args == S.One
        dimension = S.One
    return Quantity(quantity_scale[0], dimension=dimension)


dimensionless = Dimension(S.One)
