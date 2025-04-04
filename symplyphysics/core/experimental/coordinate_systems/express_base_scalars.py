from typing import TypeAlias, Mapping
from sympy import Expr, atan2, cos, sin, sqrt, Symbol as SymSymbol
from sympy.multipledispatch import dispatch

from .coordinate_systems import (
    BaseCoordinateSystem,
    CartesianCoordinateSystem,
    CylindricalCoordinateSystem,
    SphericalCoordinateSystem,
)

ScalarMapping: TypeAlias = Mapping[SymSymbol, Expr]


@dispatch(CartesianCoordinateSystem, CylindricalCoordinateSystem)
# pylint: disable-next=function-redefined
def express_base_scalars(
    old_system: CartesianCoordinateSystem,
    new_system: CylindricalCoordinateSystem,
) -> ScalarMapping:
    x, y, z = old_system.base_scalars
    rho, phi, z_ = new_system.base_scalars

    return {
        x: rho * cos(phi),
        y: rho * sin(phi),
        z: z_,
    }


@dispatch(CartesianCoordinateSystem, SphericalCoordinateSystem)  # type: ignore[no-redef]
# pylint: disable-next=function-redefined
def express_base_scalars(
    old_system: CartesianCoordinateSystem,
    new_system: SphericalCoordinateSystem,
) -> ScalarMapping:
    x, y, z = old_system.base_scalars
    r, theta, phi = new_system.base_scalars

    return {
        x: r * sin(theta) * cos(phi),
        y: r * sin(theta) * sin(phi),
        z: r * cos(theta),
    }


@dispatch(CylindricalCoordinateSystem, CartesianCoordinateSystem)  # type: ignore[no-redef]
# pylint: disable-next=function-redefined
def express_base_scalars(
    old_system: CylindricalCoordinateSystem,
    new_system: CartesianCoordinateSystem,
) -> ScalarMapping:
    rho, phi, z = old_system.base_scalars
    x, y, z_ = new_system.base_scalars

    return {
        rho: sqrt(x**2 + y**2),
        phi: atan2(y, x),
        z: z_,
    }


@dispatch(CylindricalCoordinateSystem, SphericalCoordinateSystem)  # type: ignore[no-redef]
# pylint: disable-next=function-redefined
def express_base_scalars(
    old_system: CylindricalCoordinateSystem,
    new_system: SphericalCoordinateSystem,
) -> ScalarMapping:
    rho, phi, z = old_system.base_scalars
    r, theta, phi_ = new_system.base_scalars

    return {
        rho: r * sin(theta),
        phi: phi_,
        z: r * cos(theta),
    }


@dispatch(SphericalCoordinateSystem, CartesianCoordinateSystem)  # type: ignore[no-redef]
# pylint: disable-next=function-redefined
def express_base_scalars(
    old_system: SphericalCoordinateSystem,
    new_system: CartesianCoordinateSystem,
) -> ScalarMapping:
    r, theta, phi = old_system.base_scalars
    x, y, z = new_system.base_scalars

    return {
        r: sqrt(x**2 + y**2 + z**2),
        theta: atan2(sqrt(x**2 + y**2), z),
        phi: atan2(y, x),
    }


@dispatch(SphericalCoordinateSystem, CylindricalCoordinateSystem)  # type: ignore[no-redef]
# pylint: disable-next=function-redefined
def express_base_scalars(
    old_system: SphericalCoordinateSystem,
    new_system: CylindricalCoordinateSystem,
) -> ScalarMapping:
    r, theta, phi = old_system.base_scalars
    rho, phi_, z = new_system.base_scalars

    return {
        r: sqrt(rho**2 + z**2),
        theta: atan2(rho, z),
        phi: phi_,
    }


# A fall-through case
# 1) when the type of either of the inputs hasn't been registered in `dispatch`, this will raise a `TypeError`
# 2) or when both inputs belong to the same coordinate system type, a trivial substitution is returned.
@dispatch(BaseCoordinateSystem, BaseCoordinateSystem)  # type: ignore[no-redef]
# pylint: disable-next=function-redefined
def express_base_scalars(
    old_system: BaseCoordinateSystem,
    new_system: BaseCoordinateSystem,
) -> ScalarMapping:
    old_type = type(old_system)
    new_type = type(new_system)
    if old_type is not new_type:
        raise TypeError(
            f"Conversion between {old_type.__name__} and {new_type.__name__} is not supported.")

    o1, o2, o3 = old_system.base_scalars
    n1, n2, n3 = new_system.base_scalars

    return {
        o1: n1,
        o2: n2,
        o3: n3,
    }


__all__ = ["express_base_scalars"]
