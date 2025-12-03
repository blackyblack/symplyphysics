from __future__ import annotations

from typing import Iterable

from sympy import Basic, Symbol as SymSymbol, Expr, Integral, diff

from ..vectors import VectorCross, VectorNorm
from .point import AppliedPoint
from .vector import CoordinateVector


class Surface(Basic):
    """
    Topologically, a **surface** is a manifold of dimension two. Usually we deal with **parametric
    differentiable surfaces**, which are defined by a continuous function of two variables whose
    domain is an open subset of the Euclidean plane.
    """

    @property
    def parameters(self) -> tuple[SymSymbol, SymSymbol]:
        """The symbols that parametrize the surface."""

        return self.args[0]

    @property
    def parametrization(self) -> AppliedPoint:
        """Representation of the surface as a function of the parameters."""

        return self.args[1]

    def __new__(
        cls,
        parameters: Iterable[SymSymbol] | tuple[SymSymbol, SymSymbol],
        parametrization: AppliedPoint,
    ) -> Surface:
        p1, p2 = parameters

        return super().__new__(cls, (p1, p2), parametrization)  # pylint: disable=too-many-function-args

    def to_integral(
        self,
        expr: Expr,
        first_bounds: tuple[Expr, Expr],
        second_bounds: tuple[Expr, Expr],
    ) -> Integral:
        p1, p2 = self.parameters

        lo1, hi1 = first_bounds
        lo2, hi2 = second_bounds

        return Integral(expr, (p1, lo1, hi1), (p2, lo2, hi2))

    @property
    def tangent_vectors(self) -> tuple[CoordinateVector, CoordinateVector]:
        r"""
        For any variable :math:`s`, the tangent vector along the isoline of :math:`s` is given by

        .. math::

            \vec{T}_s
            = \frac{\partial \vec{r}}{\partial s}
            = \sum_i \frac{\partial \vec{r}}{\partial q^i} \frac{\partial q^i}{\partial s}
            = \sum_i \frac{\partial q^i}{\partial s} \vec{e}_i
            = \sum_i h_i \frac{\partial q^i}{\partial s} \hat{e}_i

        where we used the chain rule, the relation :math:`\vec{e}_i = h_i \hat{e}_i`, and the
        following notation:

        * :math:`\vec{r}` is the position vector,
        * :math:`\q^i` is the :math:`i`-th curvilinear coordinate,
        * :math:`\vec{e}_i` is the :math:`i`-th covariant basis vector,
        * :math:`h_i = \left| \vec{e}_i \right|` is the :math:`i`-th Lamé coefficient,
        * :math:`\hat{e}_i` is the :math:`i`-th unit basis vector.
        """

        u, v = self.parameters
        system = self.parametrization.system
        parametrization_dict = self.parametrization.coordinates

        hs = system.lame_coefficients()
        hs = [h.subs(parametrization_dict) for h in hs]

        parametrization = parametrization_dict.values()

        t_u_components = [diff(q, u) * h for q, h in zip(parametrization, hs)]
        t_u = CoordinateVector(t_u_components, system, self.parametrization)

        t_v_components = [diff(q, v) * h for q, h in zip(parametrization, hs)]
        t_v = CoordinateVector(t_v_components, system, self.parametrization)

        return t_u, t_v

    @property
    def vector_area_multiple(self) -> CoordinateVector:
        r"""
        Returns the cross product between the tangent vectors at any point of the given surface.
        The infinitesimal vector area can be found by multiplying the result with the infinitesimal
        changes of the surface's parameters:

        .. math::

            d \vec{A} = \vec{T}_u \times \vec{T}_v \, du \, dv
        """

        t_u, t_v = self.tangent_vectors
        return VectorCross(t_u, t_v)

    @property
    def scalar_area_multiple(self) -> Expr:
        r"""
        Returns the norm of the cross product between the tangent vectors at any point of the given
        surface. The infinitesimal area of the parallelogram spanned by the tangent vectors can be
        found by multiplying the result with the infinitesimal changes of the surface's parameters:

        .. math::

            dA = \left| \vec{T}_u \times \vec{T}_v \right| du \, dv
        """

        return VectorNorm(self.vector_area_multiple)
