from typing import TypeAlias
from sympy import Expr
from symplyphysics.core.dimensions import ScalarValue

ParameterType: TypeAlias = Expr
ParameterLimits: TypeAlias = tuple[ParameterType, ScalarValue, ScalarValue]
