import functools
import inspect
import numbers
from typing import List
from sympy import S
from sympy.physics.units import Quantity, Dimension
from sympy.physics.units.systems.si import SI

from symplyphysics.core.symbols.symbols import Function, Symbol
from ..core.vectors.vectors import Vector
from .errors import UnitsError


def assert_equivalent_dimension(
        arg: Quantity, decorator_name, param_name, func_name, expected_unit: Dimension):
    if not isinstance(expected_unit, Dimension):
        raise TypeError(f"Argument '{expected_unit.name}' to decorator '{decorator_name}'"
            f" should be sympy.physics.units.Dimension.")
    #HACK: this allows to treat angle type as dimensionless
    expected_dimension = expected_unit.subs("angle", S.One)
    if isinstance(arg,
        numbers.Number) and SI.get_dimension_system().is_dimensionless(expected_dimension):
        return
    if not isinstance(arg, Quantity):
        raise TypeError(f"Argument '{param_name}' to function '{func_name}'"
            f" should be sympy.physics.units.Quantity.")

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
def validate_input(**decorator_kwargs):

    def validate_func(func):

        @functools.wraps(func)
        def wrapper_validate(*args, **kwargs):
            wrapped_signature = inspect.signature(func)
            bound_args = wrapped_signature.bind(*args, **kwargs)
            for param in wrapped_signature.parameters.values():
                if param.name in decorator_kwargs:
                    arg = bound_args.arguments[param.name]

                    components = None
                    if isinstance(arg, Vector):
                        components = arg.components
                    elif isinstance(arg, List):
                        components = arg

                    expected_unit = decorator_kwargs[param.name]
                    if isinstance(expected_unit, Symbol) or isinstance(expected_unit, Function):
                        expected_unit = expected_unit.dimension
                    if components is None:
                        assert_equivalent_dimension(arg, "validate_input", param.name,
                            func.__name__, expected_unit)
                    else:
                        for idx, c in enumerate(components):
                            assert_equivalent_dimension(c, "validate_input", f"{param.name}[{idx}]",
                                func.__name__, expected_unit)
            return func(*args, **kwargs)

        return wrapper_validate

    return validate_func


def _assert_expected_unit(value: Quantity | Vector | List,
    expected_unit: Dimension | Symbol | Function, function_name: str, validator_name: str):
    components = None
    if isinstance(value, Vector):
        components = value.components
    elif isinstance(value, List):
        components = value

    expected_dimension = expected_unit
    if isinstance(expected_unit, Symbol) or isinstance(expected_unit, Function):
        expected_dimension = expected_unit.dimension

    if components is None:
        assert_equivalent_dimension(value, validator_name, "return", function_name,
            expected_dimension)
    else:
        for idx, c in enumerate(components):
            assert_equivalent_dimension(c, validator_name, f"return[{idx}]", function_name,
                expected_dimension)


# Validates the output quantity. Output should be sympy.physics.units.Quantity, list of Quantity or
# Vector of Quantity type.
# Input should be should be Symbol with dimension property, Quantity or Dimension.
# Example:
# @validate_output(units.length**2)
# @validate_output(body_volume)
def validate_output(expected_unit):

    def validate_func(func):

        @functools.wraps(func)
        def wrapper_validate(*args, **kwargs):
            ret = func(*args, **kwargs)
            _assert_expected_unit(ret, expected_unit, func.__name__, "validate_output")
            return ret

        return wrapper_validate

    return validate_func


# Validates that output quantity has the same dimension as input quantity. Output and input parameter should be
# sympy.physics.units.Quantity, list of Quantity or Vector of Quantity type.
# Example:
# @validate_output_same("param1")
def validate_output_same(param_name):

    def validate_func(func):

        @functools.wraps(func)
        def wrapper_validate(*args, **kwargs):
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

            components = None
            if isinstance(ret, Vector):
                components = ret.components
            elif isinstance(ret, List):
                components = ret

            expected_dimension = None
            if components is None:
                expected_dimension = expected_unit.dimension
            else:
                if len(components) == 0:
                    raise UnitsError(
                        f"Argument '{param_name}' to decorator 'validate_output_same' is a "
                        f"List but does not have any components to derived expected dimension")
                expected_dimension = components[0].dimension
                for c in components:
                    if SI.get_dimension_system().equivalent_dims(c, expected_dimension):
                        continue
                    raise UnitsError(
                        f"Argument '{param_name}' to function '{func.__name__}' must"
                        f" have all component dimensions equivalent to '{expected_dimension.name}'")

            _assert_expected_unit(ret, expected_dimension, func.__name__, "validate_output_same")
            return ret

        return wrapper_validate

    return validate_func
