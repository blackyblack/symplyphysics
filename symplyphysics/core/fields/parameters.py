from typing import TypeAlias
from sympy import Expr
from ..dimensions import ScalarValue

ParameterType: TypeAlias = Expr
ParameterLimits: TypeAlias = tuple[ParameterType, ScalarValue, ScalarValue]
