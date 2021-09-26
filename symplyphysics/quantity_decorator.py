import functools
import inspect
from sympy.physics.units import Quantity, Dimension
from sympy.physics.units.systems.si import SI
from .errors import UnitsError

def assert_equivalent_dimension(arg: Quantity, decorator_name, param_name, func_name, expected_unit: Dimension):
    if not isinstance(arg, Quantity):
        raise TypeError(f"Argument '{param_name}' to function '{func_name}'"
            f" should be sympy.physics.units.Quantity.")
    if not isinstance(expected_unit, Dimension):
        raise TypeError(f"Argument '{expected_unit.name}' to decorator '{decorator_name}'"
            f" should be sympy.physics.units.Dimension.")
    if not SI.get_dimension_system().equivalent_dims(arg.dimension, expected_unit):
        raise UnitsError(f"Argument '{param_name}' to function '{func_name}' must "
            f"be in units equivalent to '{expected_unit.name}'")

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

def validate_output(expected_unit):
    def validate_func(func):
        @functools.wraps(func)
        def wrapper_validate(*args, **kwargs):
            ret = func(*args, **kwargs)
            assert_equivalent_dimension(ret, 'validate_output', 'return', func.__name__, expected_unit)
            return ret
        return wrapper_validate
    return validate_func