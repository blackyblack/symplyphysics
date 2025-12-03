from typing import TypeAlias, Any
from sympy import Expr

ParameterType: TypeAlias = Expr
ParameterLimits: TypeAlias = tuple[ParameterType, Any, Any]
