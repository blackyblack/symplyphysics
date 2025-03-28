from typing import TypeAlias, Mapping, Sequence, Any
from sympy import sqrt, sin, cos
from sympy.multipledispatch import dispatch

from ..vectors import VectorExpr
from . import (
    BaseCoordinateSystem,
    CartesianCoordinateSystem,
    CylindricalCoordinateSystem,
    SphericalCoordinateSystem,
)

VectorMapping: TypeAlias = Mapping[VectorExpr, VectorExpr]


@dispatch(CartesianCoordinateSystem, CylindricalCoordinateSystem)
# pylint: disable-next=function-redefined
def express_base_vectors(
        old_system: CartesianCoordinateSystem,
        new_system: CylindricalCoordinateSystem,
        *,
        old_args: Sequence[Any] = (),
        new_args: Sequence[Any] = (),
) -> VectorMapping:
    i, j, k = old_system.base_vectors(*old_args)
    e_rho, e_phi, e_z = new_system.base_vectors(*new_args)

    phi = new_system.phi

    return {
        i: e_rho * cos(phi) - e_phi * sin(phi),
        j: e_rho * sin(phi) + e_phi * cos(phi),
        k: e_z,
    }


@dispatch(CartesianCoordinateSystem, SphericalCoordinateSystem)  # type: ignore[no-redef]
# pylint: disable-next=function-redefined
def express_base_vectors(
        old_system: CartesianCoordinateSystem,
        new_system: SphericalCoordinateSystem,
        *,
        old_args: Sequence[Any] = (),
        new_args: Sequence[Any] = (),
) -> VectorMapping:
    i, j, k = old_system.base_vectors(*old_args)
    e_r, e_theta, e_phi = new_system.base_vectors(*new_args)

    theta = new_system.theta
    phi = new_system.phi

    return {
        i: (e_r * sin(theta) + e_theta * cos(theta)) * cos(phi) - e_phi * sin(phi),
        j: (e_r * sin(theta) + e_theta * cos(theta)) * sin(phi) + e_phi * cos(phi),
        k: e_r * cos(theta) - e_theta * sin(theta)
    }


@dispatch(CylindricalCoordinateSystem, CartesianCoordinateSystem)  # type: ignore[no-redef]
# pylint: disable-next=function-redefined
def express_base_vectors(
        old_system: CylindricalCoordinateSystem,
        new_system: CartesianCoordinateSystem,
        *,
        old_args: Sequence[Any] = (),
        new_args: Sequence[Any] = (),
) -> VectorMapping:
    e_rho, e_phi, e_z = old_system.base_vectors(*old_args)
    i, j, k = new_system.base_vectors(*new_args)

    x, y, _ = new_system.base_scalars

    rho = sqrt(x**2 + y**2)

    return {
        e_rho: (i * x + j * y) / rho,
        e_phi: (-i * y + j * x) / rho,
        e_z: k,
    }


@dispatch(CylindricalCoordinateSystem, SphericalCoordinateSystem)  # type: ignore[no-redef]
# pylint: disable-next=function-redefined
def express_base_vectors(
        old_system: CylindricalCoordinateSystem,
        new_system: SphericalCoordinateSystem,
        *,
        old_args: Sequence[Any] = (),
        new_args: Sequence[Any] = (),
) -> VectorMapping:
    e_rho, e_phi, e_z = old_system.base_vectors(*old_args)
    e_r, e_theta, e_phi_ = new_system.base_vectors(*new_args)

    theta = new_system.theta

    return {
        e_rho: e_r * sin(theta) + e_theta * cos(theta),
        e_phi: e_phi_,
        e_z: e_r * cos(theta) - e_theta * sin(theta),
    }


@dispatch(SphericalCoordinateSystem, CartesianCoordinateSystem)  # type: ignore[no-redef]
# pylint: disable-next=function-redefined
def express_base_vectors(
        old_system: SphericalCoordinateSystem,
        new_system: CartesianCoordinateSystem,
        *,
        old_args: Sequence[Any] = (),
        new_args: Sequence[Any] = (),
) -> VectorMapping:
    e_r, e_theta, e_phi = old_system.base_vectors(*old_args)
    i, j, k = new_system.base_vectors(*new_args)

    x, y, z = new_system.base_scalars

    r = sqrt(x**2 + y**2 + z**2)
    rho = sqrt(x**2 + y**2)

    return {
        e_r: (i * x + j * y + k * z) / r,
        e_theta: i * (x * z / (rho * r)) + j * (y * z / (rho * r)) - k * (rho / r),
        e_phi: (-i * y + j * x) / rho,
    }


@dispatch(SphericalCoordinateSystem, CylindricalCoordinateSystem)  # type: ignore[no-redef]
# pylint: disable-next=function-redefined
def express_base_vectors(
        old_system: SphericalCoordinateSystem,
        new_system: CylindricalCoordinateSystem,
        *,
        old_args: Sequence[Any] = (),
        new_args: Sequence[Any] = (),
) -> VectorMapping:
    e_r, e_theta, e_phi = old_system.base_vectors(*old_args)
    e_rho, e_phi_, e_z = new_system.base_vectors(*new_args)

    rho = new_system.rho
    z = new_system.z

    r = sqrt(rho**2 + z**2)

    return {
        e_r: (e_rho * rho + e_z * z) / r,
        e_theta: (e_rho * z - e_z * rho) / r,
        e_phi: e_phi_,
    }


@dispatch(BaseCoordinateSystem, BaseCoordinateSystem)  # type: ignore[no-redef]
# pylint: disable-next=function-redefined
def express_base_vectors(
        old_system: BaseCoordinateSystem,
        new_system: BaseCoordinateSystem,
        *,
        old_args: Sequence[Any] = (),
        new_args: Sequence[Any] = (),
) -> VectorMapping:
    old_type = type(old_system)
    new_type = type(new_system)
    if old_type is not new_type:
        raise TypeError(
            f"Conversion between {old_type.__name__} and {new_type.__name__} is not supported.")

    if old_system is new_system:
        return {}

    o1, o2, o3 = old_system.base_vectors(*old_args)
    n1, n2, n3 = new_system.base_vectors(*new_args)

    return {
        o1: n1,
        o2: n2,
        o3: n3,
    }


__all__ = ["express_base_vectors"]
