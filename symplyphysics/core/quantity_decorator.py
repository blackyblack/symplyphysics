import functools
import inspect
from sympy.core.singleton import S
from sympy.physics.units import Quantity, Dimension
from sympy.physics.units.systems.si import SI
from .errors import UnitsError

def assert_equivalent_dimension(arg: Quantity, decorator_name, param_name, func_name, expected_unit: Dimension):
    if not isinstance(arg, Quantity):
        raise TypeError(f"Argument '{param_name}' to function '{func_name}'"
            f" should be sympy.physics.units.Quantity.")

    scale_factor = SI.get_quantity_scale_factor(arg)
    # zero can be of any dimension
    if scale_factor == S.Zero:
        return

    if not isinstance(expected_unit, Dimension):
        raise TypeError(f"Argument '{expected_unit.name}' to decorator '{decorator_name}'"
            f" should be sympy.physics.units.Dimension.")
    if not SI.get_dimension_system().equivalent_dims(arg.dimension, expected_unit):
        raise UnitsError(f"Argument '{param_name}' to function '{func_name}' must "
            f"be in units equivalent to '{expected_unit.name}'")
    if scale_factor.free_symbols:
        raise UnitsError(f"Argument '{param_name}' to function '{func_name}' should "
            f"not contain free symbols")

# Validates the input quantities. Input parameters should be sympy.physics.units.Quantity type.
# Example:
# @validate_input(param1_=units.length, param2_=(1 / units.length))
def validate_input(**decorator_kwargs):
    def validate_func(func):
        @functools.wraps(func)
        def wrapper_validate(*args, **kwargs):
            wrapped_signature = inspect.signature(func)
            bound_args = wrapped_signature.bind(*args, **kwargs)
            for param in wrapped_signature.parameters.values():
                if param.name in decorator_kwargs:
                    expected_unit = decorator_kwargs[param.name]
                    arg = bound_args.arguments[param.name]
                    assert_equivalent_dimension(arg, 'validate_input', param.name, func.__name__, expected_unit)
            return func(*args, **kwargs)
        return wrapper_validate
    return validate_func

# Validates the output quantity. Output should be sympy.physics.units.Quantity type.
# Example:
# @validate_output(units.length**2)
def validate_output(expected_unit):
    def validate_func(func):
        @functools.wraps(func)
        def wrapper_validate(*args, **kwargs):
            ret = func(*args, **kwargs)
            assert_equivalent_dimension(ret, 'validate_output', 'return', func.__name__, expected_unit)
            return ret
        return wrapper_validate
    return validate_func

# Validates that output quantity has the same dimension as input quantity. Output and input parameter should be sympy.physics.units.Quantity type.
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
            assert_equivalent_dimension(ret, 'validate_output_same', param.name, func.__name__, expected_unit.dimension)
            return ret
        return wrapper_validate
    return validate_func