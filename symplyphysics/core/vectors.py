from typing import List
from sympy.vector import CoordSys3D, Vector
from sympy.vector.operators import _get_coord_systems


# Converts Python List to SymPy Vector
def array_to_sympy_vector(coord_system_: CoordSys3D, array_: List) -> Vector:
    result_vector = Vector.zero
    base_vectors = coord_system_.base_vectors()
    for idx in range(min(len(base_vectors), len(array_))):
        result_vector = result_vector + base_vectors[idx] * array_[idx]
    return result_vector

# Converts SymPy Vector to Python List
# SymPy vector is an expression that looks like C.i + C.j, where C is CoordSys3D
def sympy_vector_to_array(sympy_vector_: Vector) -> List:
    result_vector = []
    coord_system_set = _get_coord_systems(sympy_vector_)
    if len(coord_system_set) == 0: return []
    if len(coord_system_set) > 1:
        coord_sys_names = [str(c) for c in coord_system_set]
        raise TypeError(f"Different coordinate systems in expression: {str(coord_sys_names)}")
    coord_system = next(iter(coord_system_set))
    if not sympy_vector_.is_Vector:
        raise TypeError(f"Expression should be Vector")
    as_matrix = sympy_vector_.to_matrix(coord_system)
    for e in as_matrix:
        result_vector.append(e)
    return result_vector
