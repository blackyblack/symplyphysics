from typing import TypeAlias, Callable
from sympy import Expr

WaveFunction: TypeAlias = Callable[[Expr, Expr], Expr]
Observable: TypeAlias = Callable[[WaveFunction], WaveFunction]
