from __future__ import annotations
from typing import Optional, Sequence
from sympy import Expr
from sympy.vector import express, Vector as SymVector
from sympy.vector.operators import _get_coord_systems
from sympy.physics.units import Dimension
from sympy.physics.units.definitions.dimension_definitions import angle as angle_type

from ..dimensions import assert_equivalent_dimension, dimensionless, ScalarValue
from ..symbols.quantities import Quantity, subs_list
from ..symbols.symbols import DimensionSymbol, next_name
from ..coordinate_systems.coordinate_systems import CoordinateSystem


# Contains list of SymPy expressions or any numbers as components.
# Contains coordinate system to prevent using vector arithmetics with non compatible
# coordinate systems.
# Vector assumes to have [0, 0, 0] point of origin. All vector related laws work with this assumption
# add do not differentiate between position shifted vectors. But it does not allow to properly rebase
# vectors to another coordinate system.
# Therefore for all physical applications vectors should assume various origin point and should be
# defined dynamically, eg [C.x, C.y] or [parameter1, parameter2].
#NOTE: use Python lists when coordinate system does not matter
class Vector:
    #NOTE: 4 and higher dimensional vectors are not supported cause of using CoordSys3D
    #      to allow rebasing vector coordinate system.
    _coordinate_system: CoordinateSystem
    _components: list[ScalarValue]

    def __init__(
        self,
        components: Sequence[ScalarValue],
        coordinate_system: CoordinateSystem = CoordinateSystem(CoordinateSystem.System.CARTESIAN)
    ) -> None:
        self._coordinate_system = coordinate_system
        self._components = list(components)

    @property
    def coordinate_system(self) -> CoordinateSystem:
        return self._coordinate_system

    @property
    def components(self) -> Sequence[ScalarValue]:
        return self._components

    # Converts SymPy Vector to Vector
    # SymPy vector is an expression that looks like C.i + C.j, where C is CoordSys3D
    @staticmethod
    def from_sympy_vector(sympy_vector_: SymVector, coordinate_system: CoordinateSystem) -> Vector:
        if sympy_vector_ == SymVector.zero:
            return Vector([], coordinate_system)
        coord_system_set = _get_coord_systems(sympy_vector_)
        coord_system = None
        if len(coord_system_set) > 1:
            coord_sys_names = [str(c) for c in coord_system_set]
            raise TypeError(f"Different coordinate systems in expression: {str(coord_sys_names)}")
        if len(coord_system_set) > 0:
            coord_system = next(iter(coord_system_set))
            if coord_system != coordinate_system.coord_system:
                raise TypeError(
                    f"Different coordinate systems in expression and argument: {str(coord_system)} vs {str(coordinate_system.coord_system)}"
                )
        components = list(sympy_vector_.to_matrix(coord_system))
        return Vector(components, coordinate_system)

    # Converts Vector to SymPy Vector
    def to_sympy_vector(self) -> SymVector:
        result_vector = SymVector.zero
        base_vectors = self.coordinate_system.coord_system.base_vectors()
        for idx in range(min(len(base_vectors), len(self.components))):
            result_vector = result_vector + base_vectors[idx] * self.components[idx]
        return result_vector

    # Convert vector coordinate system to new basis and construct new vector.
    # Rebased vector should be the same as old vector but in new coordinate system.
    def rebase(self, coordinate_system: CoordinateSystem) -> Vector:
        vector_: Vector = self
        if self.coordinate_system.coord_system_type != coordinate_system.coord_system_type:
            new_scalars = list(
                self.coordinate_system.transformation_to_system(
                coordinate_system.coord_system_type))
            # now take each component of vector and assign them to base_scalars, eg x, y, z
            # replace each component of new_scalars with assigned x, y, z, eg r -> x, theta -> y
            # build new vector from these components
            for i, scalar in enumerate(self.coordinate_system.coord_system.base_scalars()):
                new_component = 0 if i >= len(self.components) else self.components[i]
                for j, old_scalar in enumerate(new_scalars):
                    new_scalars[j] = old_scalar.subs(scalar, new_component)
            vector_ = Vector(new_scalars, self.coordinate_system)
        # We do not want to maintain own vector transformation functions, so
        # we convert our vector to SymPy format, transform it and convert back to Vector.
        sympy_vector = vector_.to_sympy_vector()
        transformed_vector_sympy = express(sympy_vector,
            coordinate_system.coord_system,
            None,
            variables=True)
        return Vector.from_sympy_vector(transformed_vector_sympy, coordinate_system)


class QuantityVector(DimensionSymbol):
    # Vector of scalar values, ie expressions or numbers
    _inner_vector: Vector

    def __init__(self,
        components: Sequence[Quantity | ScalarValue],
        coordinate_system: CoordinateSystem = CoordinateSystem(CoordinateSystem.System.CARTESIAN),
        *,
        dimension: Optional[Dimension] = None) -> None:
        quantities = [
            c if isinstance(c, Quantity) else Quantity(c, dimension=dimension) for c in components
        ]
        # find first dimension with non-zero scale factor
        if dimension is None:
            dimension = dimensionless
            for q in quantities:
                if q.scale_factor != 0:
                    dimension = q.dimension
                    break
        scale_factors = []
        for idx, c in enumerate(quantities):
            dimension_to_check = angle_type if CoordinateSystem.is_angle_component(
                coordinate_system.coord_system_type, idx) else dimension
            assert_equivalent_dimension(c, c.display_name, "QuantityVector", dimension_to_check)
            scale_factors.append(c.scale_factor)
        DimensionSymbol.__init__(self, next_name("VEC"), dimension)
        self._inner_vector = Vector(scale_factors, coordinate_system)

    @property
    def components(self) -> Sequence[Quantity]:
        quantities = []
        for idx, c in enumerate(self._inner_vector.components):
            dimension_to_set = angle_type if CoordinateSystem.is_angle_component(
                self._inner_vector.coordinate_system.coord_system_type, idx) else self.dimension
            quantities.append(Quantity(c, dimension=dimension_to_set))
        return quantities

    @property
    def coordinate_system(self) -> CoordinateSystem:
        return self._inner_vector.coordinate_system

    def to_base_vector(self) -> Vector:
        # Quantities are expressions - do not need additional sympify
        return Vector(self.components, self.coordinate_system)

    @staticmethod
    def from_base_vector(vector: Vector, *,
        dimension: Optional[Dimension] = None,
        subs: Optional[dict[Expr, Quantity]] = None) -> QuantityVector:
        components = vector.components
        if subs is not None:
            components = subs_list(components, subs)
        return QuantityVector(components, vector.coordinate_system, dimension=dimension)
