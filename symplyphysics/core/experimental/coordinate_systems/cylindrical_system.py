from typing import Optional, Sequence

from sympy import Expr, S, cos, sin, sqrt, atan2, Q
from sympy.logic.boolalg import Boolean

from symplyphysics import Symbol, units, angle_type

from .base_system import BaseCoordinateSystem


class CylindricalCoordinateSystem(BaseCoordinateSystem):

    def generate_base_scalars(self) -> tuple[Symbol, Symbol, Symbol]:
        return (
            Symbol("rho", units.length, display_latex="\\rho", nonnegative=True),
            Symbol("phi", angle_type, display_latex="\\varphi", real=True),
            Symbol("z", units.length, real=True),
        )

    def cartesian_transform(self, base_scalars: Optional[Sequence[Expr]] = None) -> Sequence[Expr]:
        rho, phi, z = base_scalars or self.base_scalars

        x = rho * cos(phi)
        y = rho * sin(phi)

        return x, y, z

    def inverse_cartesian_transform(self, cartesian_scalars: Sequence[Expr]) -> Sequence[Expr]:
        x, y, z = cartesian_scalars

        rho = sqrt(x**2 + y**2)
        phi = atan2(y, x)

        return rho, phi, z

    def assumption(self, base_scalars: Optional[Sequence[Expr]] = None) -> Boolean:
        rho = (base_scalars or self.base_scalars)[0]

        return Q.positive(rho)  # pylint: disable=too-many-function-args

    def generate_lame_coefficients(self) -> tuple[Expr, Expr, Expr]:
        rho, _, _ = self.base_scalars

        h_rho = S.One
        h_phi = rho
        h_z = S.One

        return h_rho, h_phi, h_z


__all__ = ["CylindricalCoordinateSystem"]
