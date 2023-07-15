import functools
import inspect
from typing import Any, Callable, Sequence
from sympy import S
from sympy.physics.units import Quantity as SymQuantity, Dimension
from sympy.physics.units.systems.si import SI

from symplyphysics.core.symbols.quantities import Quantity

from ..core.symbols.symbols import DimensionSymbol, Function, Symbol
from ..core.vectors.vectors import Vector
from .errors import UnitsError


def assert_equivalent_dimension(
        arg: SymQuantity | float, param_name: str, func_name: str, expected_unit: Dimension):
    #HACK: this allows to treat angle type as dimensionless
    expected_dimension = expected_unit.subs("angle", S.One)
    if isinstance(arg, (float | int)):
        if SI.get_dimension_system().is_dimensionless(expected_dimension):
            return
        raise TypeError(f"Argument '{param_name}' to function '{func_name}'"
            f" is Number but '{expected_dimension}' is not dimensionless")
    arg_quantity = Quantity(arg)
    # zero can be of any dimension
    if arg_quantity.scale_factor == S.Zero:
        return
    #HACK: this allows to treat angle type as dimensionless
    arg_dimension = arg_quantity.dimension.subs("angle", S.One)
    # angle is dimensionless but equivalent_dims() fails to compare it
    if SI.get_dimension_system().is_dimensionless(
            expected_dimension) and SI.get_dimension_system().is_dimensionless(arg_dimension):
        return
    if not SI.get_dimension_system().equivalent_dims(arg_dimension, expected_dimension):
        raise UnitsError(f"Argument '{param_name}' to function '{func_name}' must "
            f"be in units equivalent to '{expected_dimension.name}'")
    if arg_quantity.scale_factor.free_symbols:
        raise UnitsError(f"Argument '{param_name}' to function '{func_name}' should "
            f"not contain free symbols")


def _assert_expected_unit(value: SymQuantity | Vector | Sequence[SymQuantity | float],
    expected_units: Dimension | Symbol | Function | Sequence[Dimension | Symbol | Function],
    param_name: str, function_name: str):
    components: list[SymQuantity | float] = []
    indexed = False
    if isinstance(value, Vector):
        for c in value.components:
            if not isinstance(c, (SymQuantity | float | int)):
                raise TypeError(f"Argument '{param_name}' to function '{function_name}'"
                    f" is Vector but its component {c} is not Quantity or Number")
            components.append(c)
        indexed = True
    elif isinstance(value, Sequence):
        components = list(value)
        indexed = True
    else:
        components.append(value)

    is_tuple = False
    if isinstance(expected_units, Sequence):
        is_tuple = True
    else:
        expected_units = [expected_units]

    expected_unit_dimensions: list[Dimension] = []
    for u in list(expected_units):
        d = u.dimension if isinstance(u, DimensionSymbol) else u
        expected_unit_dimensions.append(d)

    for idx, c in enumerate(components):
        param_name_indexed = f"{param_name}[{idx}]" if indexed else param_name
        expected_dimension = expected_unit_dimensions[
            idx] if is_tuple else expected_unit_dimensions[0]
        assert_equivalent_dimension(c, param_name_indexed, function_name, expected_dimension)


# Validates the input quantities. Input parameters should be sympy.physics.units.Quantity, list of Quantity or
# Vector of Quantity type.
# Unit should be should be Symbol with dimension property, or Dimension.
# Example:
# @validate_input(param1_=units.length, param2_=(1 / units.length))
# @validate_input(param1_=body_mass, param2_=body_volume)
def validate_input(**decorator_kwargs: Any) -> Callable[[Callable[..., Any]], Callable[..., Any]]:

    def validate_func(func: Callable[..., Any]) -> Callable[..., Any]:

        @functools.wraps(func)
        def wrapper_validate(*args: Any, **kwargs: Any) -> Any:
            wrapped_signature = inspect.signature(func)
            bound_args = wrapped_signature.bind(*args, **kwargs)
            for param in wrapped_signature.parameters.values():
                if param.name in decorator_kwargs:
                    arg = bound_args.arguments[param.name]
                    _assert_expected_unit(arg, decorator_kwargs[param.name], param.name,
                        func.__name__)
            return func(*args, **kwargs)

        return wrapper_validate

    return validate_func


# Validates the output quantity. Output should be sympy.physics.units.Quantity, list of Quantity or
# Vector of Quantity type.
# Input should be should be Symbol with dimension property, Quantity or Dimension.
# Example:
# @validate_output(units.length)
# @validate_output(body_volume)
def validate_output(
        expected_unit: Dimension | Symbol | Function) -> Callable[[Any], Callable[..., Any]]:

    def validate_func(func: Callable[..., Any]) -> Callable[..., Any]:

        @functools.wraps(func)
        def wrapper_validate(*args: Any, **kwargs: Any) -> Any:
            ret = func(*args, **kwargs)
            _assert_expected_unit(ret, expected_unit, "return", func.__name__)
            return ret

        return wrapper_validate

    return validate_func


# Validates that output quantity has the same dimension as input quantity. Output and input parameter should be
# sympy.physics.units.Quantity, list of Quantity or Vector of Quantity type.
# Example:
# @validate_output_same("param1")
def validate_output_same(param_name: str) -> Callable[[Any], Callable[..., Any]]:

    def validate_func(func: Callable[..., Any]) -> Callable[..., Any]:

        @functools.wraps(func)
        def wrapper_validate(*args: Any, **kwargs: Any) -> Any:
            wrapped_signature = inspect.signature(func)
            bound_args = wrapped_signature.bind(*args, **kwargs)
            expected_unit = None
            for param in wrapped_signature.parameters.values():
                if param.name == param_name:
                    arg = bound_args.arguments[param.name]
                    expected_unit = arg
                    break
            if expected_unit is None:
                raise TypeError(f"Argument '{param_name}' to decorator 'validate_output_same'"
                    f" should be in function parameters")
            ret = func(*args, **kwargs)

            _assert_expected_unit(ret, expected_unit, "return", func.__name__)
            return ret

        return wrapper_validate

    return validate_func
