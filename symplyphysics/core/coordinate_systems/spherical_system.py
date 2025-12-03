from typing import Optional, Sequence, Iterable

from sympy import Expr, S, cos, sin, sqrt, atan2, Q, tan
from sympy.logic.boolalg import Boolean
from sympy.physics import units
from sympy.physics.units.definitions.dimension_definitions import angle as angle_type

from symplyphysics.core.symbols.symbols import Symbol

from .base_system import BaseCoordinateSystem


class SphericalCoordinateSystem(BaseCoordinateSystem):

    def generate_base_scalars(self) -> tuple[Symbol, Symbol, Symbol]:
        return (
            Symbol("r", units.length, nonnegative=True),
            Symbol("theta", angle_type, display_latex="\\theta", nonnegative=True),
            Symbol("phi", angle_type, display_latex="\\varphi", real=True),
        )

    def cartesian_transform(self, base_scalars: Optional[Sequence[Expr]] = None) -> Sequence[Expr]:
        r, theta, phi = base_scalars or self.base_scalars

        x = r * sin(theta) * cos(phi)
        y = r * sin(theta) * sin(phi)
        z = r * cos(theta)

        return x, y, z

    def inverse_cartesian_transform(self, cartesian_scalars: Sequence[Expr]) -> Sequence[Expr]:
        x, y, z = cartesian_scalars

        r = sqrt(x**2 + y**2 + z**2)
        theta = atan2(sqrt(x**2 + y**2), z)
        phi = atan2(y, x)

        return r, theta, phi

    def assumption(self, base_scalars: Optional[Sequence[Expr]] = None) -> Boolean:
        r, theta, *_ = base_scalars or self.base_scalars

        # pylint: disable-next=too-many-function-args
        return Q.positive(r) & Q.positive(sin(theta)) & Q.positive(tan(theta))

    def generate_lame_coefficients(self) -> tuple[Expr, Expr, Expr]:
        r, theta, _ = self.base_scalars

        h_r = S.One
        h_theta = r
        h_phi = r * sin(theta)

        return h_r, h_theta, h_phi

    def extract_position_vector_components(self, coordinates: Iterable[Expr]) -> Sequence[Expr]:
        r, _, _ = coordinates

        return r, S.Zero, S.Zero


__all__ = ["SphericalCoordinateSystem"]
