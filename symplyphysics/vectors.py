from typing import List
from sympy.vector import CoordSys3D, Vector


# Converts Python List to SymPy Vector.
def array_to_sympy_vector(coord_system_: CoordSys3D, array_: List) -> Vector:
    result_vector = Vector.zero
    array_len = len(array_)
    for idx, coord in enumerate(coord_system_):
        if idx >= array_len: break
        result_vector = result_vector + coord * array_[idx]
    return result_vector