from __future__ import annotations
from functools import partial
from typing import Any, Optional
from sympy import Expr, S, Derivative, Function
from sympy.core.add import Add
from sympy.core.mul import Mul
from sympy.core.power import Pow
from sympy.physics.units import Dimension, Quantity as SymQuantity
from sympy.physics.units.systems.si import SI
from .symbols import DimensionSymbol, next_name


class Quantity(DimensionSymbol, SymQuantity):

    def __new__(cls,
        _scale: Expr = S.One,
        *,
        display_name: Optional[str] = None,
        **assumptions: Any) -> SymQuantity:
        name = next_name("QU") if display_name is None else next_name(display_name)
        obj = SymQuantity.__new__(cls, name, None, None, None, None, None, False, **assumptions)
        return obj

    def __init__(self,
        scale: Expr = S.One,
        *,
        display_name: Optional[str] = None,
        dimension: Optional[Dimension] = None,
        **_assumptions: Any):
        (_, dimension_) = collect_factor_and_dimension(scale)
        dimension = dimension_ if dimension is None else dimension
        super().__init__(self.name if display_name is None else display_name, dimension)
        SI.set_quantity_dimension(self, dimension)
        SI.set_quantity_scale_factor(self, scale)

    # This is required for integration to work properly
    @property
    def func(self) -> partial:
        return partial(Quantity.identity, self)

    def identity(self, *_args: Any) -> Quantity:
        return self


# HACK: copy of SI._collect_factor_and_dimension with fixed exp() dimension evaluation
def collect_factor_and_dimension(expr: Expr) -> tuple[Expr, Dimension]:
    """
    Return tuple with scale factor expression and dimension expression.
    """
    if isinstance(expr, SymQuantity):
        return expr.scale_factor, expr.dimension
    elif isinstance(expr, Mul):
        factor = S.One
        dimension = Dimension(S.One)
        for arg in expr.args:
            arg_factor, arg_dim = collect_factor_and_dimension(arg)
            factor *= arg_factor
            dimension *= arg_dim
        return factor, dimension
    elif isinstance(expr, Pow):
        factor, dim = collect_factor_and_dimension(expr.base)
        exp_factor, exp_dim = collect_factor_and_dimension(expr.exp)
        if SI.get_dimension_system().is_dimensionless(exp_dim):
            exp_dim = 1
        return factor**exp_factor, dim**(exp_factor * exp_dim)
    elif isinstance(expr, Add):
        factor, dim = collect_factor_and_dimension(expr.args[0])
        for addend in expr.args[1:]:
            addend_factor, addend_dim = collect_factor_and_dimension(addend)
            # automatically convert zero to the dimension of it's additives
            if dim != addend_dim and not SI.get_dimension_system().equivalent_dims(dim, addend_dim):
                if factor == S.Zero:
                    dim = addend_dim
                elif addend_factor == S.Zero:
                    addend_dim = dim
            if dim != addend_dim and not SI.get_dimension_system().equivalent_dims(dim, addend_dim):
                raise ValueError(f"Dimension of '{addend}' is {addend_dim}, but it should be {dim}")
            factor += addend_factor
        return factor, dim
    elif isinstance(expr, Derivative):
        factor, dim = collect_factor_and_dimension(expr.args[0])
        for independent, count in expr.variable_count:
            ifactor, idim = collect_factor_and_dimension(independent)
            factor /= ifactor**count
            dim /= idim**count
        return factor, dim
    elif isinstance(expr, Function):
        fds = [collect_factor_and_dimension(arg) for arg in expr.args]
        dims = [
            Dimension(S.One) if SI.get_dimension_system().is_dimensionless(d[1]) else d[1]
            for d in fds
        ]
        return (expr.func(*(f[0] for f in fds)), *dims)
    elif isinstance(expr, Dimension):
        return S.One, expr
    else:
        return expr, Dimension(S.One)


Dimensionless = Dimension(S.One)
