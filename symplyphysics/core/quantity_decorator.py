import functools
import inspect
from typing import Any, Callable, Sequence, Optional
from sympy import S
from sympy.physics.units import Quantity as SymQuantity, Dimension
from sympy.physics.units.systems.si import SI

from ..core.symbols.symbols import DimensionSymbol, Function, Symbol
from ..core.vectors.vectors import Vector
from .errors import UnitsError


def assert_equivalent_dimension(
        arg: SymQuantity | float, param_name: str, func_name: str, expected_unit: Dimension):
    #HACK: this allows to treat angle type as dimensionless
    expected_dimension = expected_unit.subs("angle", S.One)
    if not isinstance(arg, SymQuantity):
        if SI.get_dimension_system().is_dimensionless(expected_dimension):
            return
        raise TypeError(f"Argument '{param_name}' to function '{func_name}'"
            f" is Number but '{expected_dimension}' is not dimensionless")
    scale_factor = SI.get_quantity_scale_factor(arg)
    # zero can be of any dimension
    if scale_factor == S.Zero:
        return
    #HACK: this allows to treat angle type as dimensionless
    arg_dimension = arg.dimension.subs("angle", S.One)
    # angle is dimensionless but equivalent_dims() fails to compare it
    if SI.get_dimension_system().is_dimensionless(
            expected_dimension) and SI.get_dimension_system().is_dimensionless(arg_dimension):
        return
    if not SI.get_dimension_system().equivalent_dims(arg_dimension, expected_dimension):
        raise UnitsError(f"Argument '{param_name}' to function '{func_name}' must "
            f"be in units equivalent to '{expected_dimension.name}'")
    if scale_factor.free_symbols:
        raise UnitsError(f"Argument '{param_name}' to function '{func_name}' should "
            f"not contain free symbols")


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

                    components: Optional[list[SymQuantity | float]] = None
                    if isinstance(arg, Vector):
                        components = arg.components
                    elif isinstance(arg, Sequence):
                        components = list(arg)

                    expected_unit = decorator_kwargs[param.name]
                    if isinstance(expected_unit, DimensionSymbol):
                        expected_unit = expected_unit.dimension
                    if components is None:
                        assert_equivalent_dimension(arg, param.name, func.__name__, expected_unit)
                    else:
                        for idx, c in enumerate(components):
                            assert_equivalent_dimension(c, f"{param.name}[{idx}]", func.__name__,
                                expected_unit)
            return func(*args, **kwargs)

        return wrapper_validate

    return validate_func


def _assert_expected_unit(value: SymQuantity | Vector | Sequence,
    expected_unit: Dimension | Symbol | Function, function_name: str):
    components: Optional[list[SymQuantity | float]] = None
    value_typed: Optional[SymQuantity | float] = None
    if isinstance(value, Vector):
        components = value.components
    elif isinstance(value, Sequence):
        components = list(value)
    else:
        # Make linter happy
        value_typed = value

    # Make linter happy
    expected_dimension_typed: Dimension
    if isinstance(expected_unit, DimensionSymbol):
        expected_dimension_typed = expected_unit.dimension
    else:
        expected_dimension_typed = expected_unit

    if value_typed is not None:
        assert_equivalent_dimension(value_typed, "return", function_name, expected_dimension_typed)
    if components is not None:
        for idx, c in enumerate(components):
            assert_equivalent_dimension(c, f"return[{idx}]", function_name,
                expected_dimension_typed)


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
            _assert_expected_unit(ret, expected_unit, func.__name__)
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

            components: Optional[list[SymQuantity | float]] = None
            value_typed: Optional[SymQuantity | float] = None
            if isinstance(ret, Vector):
                components = ret.components
            elif isinstance(ret, Sequence):
                components = list(ret)
            else:
                # Make linter happy
                value_typed = ret

            # Make linter happy
            expected_dimension_typed: Dimension
            if isinstance(expected_unit, DimensionSymbol):
                expected_dimension_typed = expected_unit.dimension
            else:
                expected_dimension_typed = expected_unit

            if value_typed is not None:
                _assert_expected_unit(ret, expected_dimension_typed, func.__name__)
                return ret

            assert components is not None, "Should have either components or value_typed set"

            if len(components) == 0:
                raise UnitsError(
                    f"Argument '{param_name}' to decorator 'validate_output_same' is a "
                    f"List but does not have any components to derived expected dimension")
            for c in components:
                if SI.get_dimension_system().equivalent_dims(c, expected_dimension_typed):
                    continue
                raise UnitsError(f"Argument '{param_name}' to function '{func.__name__}' must"
                    f" have all component dimensions equivalent to '{expected_dimension_typed.name}'"
                                )
            return ret

        return wrapper_validate

    return validate_func
