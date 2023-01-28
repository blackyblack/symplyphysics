from typing import Any, List
from sympy import Expr
from sympy.vector import CoordSys3D, Vector as SympyVector, express
from sympy.vector.operators import _get_coord_systems


# Detects coordinate system being used in SymPy Vector.
# Throws an exception if multiple coordinate systems are detected.
def extract_coord_system_from_sympy_vector(sympy_vector_: SympyVector) -> CoordSys3D:
    if not isinstance(sympy_vector_, Expr): return None
    coord_system_set = _get_coord_systems(sympy_vector_)
    if len(coord_system_set) > 1:
        coord_sys_names = [str(c) for c in coord_system_set]
        raise TypeError(f"Different coordinate systems in expression: {str(coord_sys_names)}")
    return None if len(coord_system_set) == 0 else next(iter(coord_system_set))

# Contains list of SymPy expressions or any numbers as components.
# Contains coordinate system to prevent using vector arithmetics with non compatible
# coordinate systems.
class Vector:
    _components: List[Any] = []
    #NOTE: 4 and higher dimensional vectors are not supported cause of using CoordSys3D
    #      to allow vector coordinate system rebasing.
    _coord_system: CoordSys3D = None

    def __init__(self, components: List[Any]=[], coord_system: CoordSys3D=None):
        self._components = components
        self._coord_system = coord_system

    @property
    def coord_system(self) -> CoordSys3D:
        return self._coord_system

    @property
    def components(self):
        return self._components

# Converts SymPy Vector to Vector
# SymPy vector is an expression that looks like C.i + C.j, where C is CoordSys3D
def vector_from_sympy_vector(sympy_vector_: Any) -> Vector:
    if sympy_vector_ == SympyVector.zero: return Vector([])
    coord_system = extract_coord_system_from_sympy_vector(sympy_vector_)
    # Do not convert if cannot detect SymPy Vector class
    if coord_system is None or not sympy_vector_.is_Vector: return Vector([sympy_vector_])
    as_matrix = sympy_vector_.to_matrix(coord_system)
    components = []
    for e in as_matrix:
        components.append(e)
    return Vector(components, coord_system)

# Converts Vector to SymPy Vector
def sympy_vector_from_vector(vector_: Vector) -> SympyVector:
    result_vector = SympyVector.zero
    if vector_.coord_system is None: return result_vector
    base_vectors = vector_.coord_system.base_vectors()
    for idx in range(min(len(base_vectors), len(vector_.components))):
        result_vector = result_vector + base_vectors[idx] * vector_.components[idx]
    return result_vector

# Convert vector coordinate system to new basis and construct new vector.
# Rebased vector should be the same as old vector but in new coordinate system.
def vector_rebase(vector_: Vector, coord_system: CoordSys3D=None) -> Vector:
    # Simply set new coordinate system if vector cannot be rebased
    if coord_system is None or vector_.coord_system is None:
        return Vector(vector_.components, coord_system)

    # We do not want to maintain own vector transformation functions, so
    # we convert our vector to SymPy format, transform it and convert back to Vector.
    sympy_vector = sympy_vector_from_vector(vector_)
    transformed_vector = express(sympy_vector, coord_system, variables=True)
    return vector_from_sympy_vector(transformed_vector)
