from pytest import raises
from sympy import S, evaluate, sqrt, Min, diff, sympify
from sympy.physics.units.systems.si import dimsys_SI
from symplyphysics import symbols, units, clone_as_symbol, clone_as_function, Quantity
from symplyphysics.core.symbols.symbols import clone_as_indexed
from symplyphysics.core.dimensions import collect_dimension


def test_mul() -> None:
    expr = symbols.position / symbols.time
    dim = collect_dimension(expr)
    assert dimsys_SI.equivalent_dims(dim, units.length / units.time)

    expr = symbols.time * symbols.temporal_frequency
    dim = collect_dimension(expr)
    assert dimsys_SI.is_dimensionless(dim)

    with evaluate(False):
        expr = S.Zero * symbols.time
    dim = collect_dimension(expr)
    assert dimsys_SI.is_dimensionless(dim)

    with evaluate(False):
        expr = symbols.time * S.Zero
    dim = collect_dimension(expr)
    assert dimsys_SI.is_dimensionless(dim)


def test_pow() -> None:
    with evaluate(False):
        expr = S.Zero**symbols.positive_number
    dim = collect_dimension(expr)
    assert dimsys_SI.is_dimensionless(dim)

    with evaluate(False):
        expr = symbols.time**S.Zero
    dim = collect_dimension(expr)
    assert dimsys_SI.is_dimensionless(dim)

    # the exponent can be a number
    expr = (symbols.mass / symbols.amount_of_substance)**2
    dim = collect_dimension(expr)
    assert dimsys_SI.equivalent_dims(dim, (units.mass / units.amount_of_substance)**2)

    # the exponent is allowed to be a single `Symbol`
    n = symbols.positive_number
    expr = symbols.current**n
    dim = collect_dimension(expr)
    assert dimsys_SI.equivalent_dims(dim, units.current**n)

    expr = symbols.current**(2 * n)
    with raises(TypeError):
        collect_dimension(expr)

    # NOTE: refer to `_collect_pow` in `collect_dimension`
    expr = (symbols.acceleration * symbols.mass)**(2 * symbols.time * symbols.temporal_frequency)
    with raises(TypeError):
        collect_dimension(expr)

    # NOTE: refer to `_collect_pow` in `collect_dimension`
    bad_expr = symbols.mass**symbols.amount_of_substance
    with raises(ValueError):
        collect_dimension(bad_expr)


def test_add() -> None:
    first_mass = clone_as_symbol(symbols.mass, subscript="1")
    second_mass = clone_as_symbol(symbols.mass, subscript="2")
    third_mass = clone_as_symbol(symbols.mass, subscript="3")

    expr = first_mass + second_mass + third_mass
    dim = collect_dimension(expr)
    assert dimsys_SI.equivalent_dims(dim, units.mass)

    expr = first_mass + second_mass + symbols.time
    with raises(ValueError):
        collect_dimension(expr)


def test_abs() -> None:
    expr = abs(symbols.time)
    dim = collect_dimension(expr)
    assert dimsys_SI.equivalent_dims(dim, units.time)

    expr = abs(sqrt(symbols.area))
    dim = collect_dimension(expr)
    assert dimsys_SI.equivalent_dims(dim, units.length)

    # the error propagates from the argument
    expr = abs(symbols.time + symbols.length)
    with raises(ValueError):
        collect_dimension(expr)


def test_min() -> None:
    first_mass = clone_as_symbol(symbols.mass, subscript="1")
    second_mass = clone_as_symbol(symbols.mass, subscript="2")
    third_mass = clone_as_symbol(symbols.mass, subscript="3")

    expr = Min(first_mass, second_mass, third_mass)
    dim = collect_dimension(expr)
    assert dimsys_SI.equivalent_dims(dim, units.mass)

    expr = Min(first_mass, second_mass, symbols.time)
    with raises(ValueError):
        collect_dimension(expr)


def test_derivative() -> None:
    position = symbols.position
    time = symbols.time
    current = clone_as_function(symbols.current, [position, time])

    expr = diff(current(position, time), position, (time, 2))
    dim = collect_dimension(expr)
    assert dimsys_SI.equivalent_dims(dim, units.current / (units.length * units.time**2))

    expr = diff(current(position, time), position)
    dim = collect_dimension(expr)
    assert dimsys_SI.equivalent_dims(dim, units.current / units.length)


def test_has_dimension_field() -> None:
    # sympy.Quantity
    expr = Quantity(4.0 * units.meter)
    dim = collect_dimension(expr)
    assert dimsys_SI.equivalent_dims(dim, units.length)

    # SymbolNew
    dim = collect_dimension(symbols.power)
    assert dimsys_SI.equivalent_dims(dim, units.power)

    # FunctionNew
    expr = clone_as_function(symbols.density, [symbols.position, symbols.time])
    dim = collect_dimension(expr)
    assert dimsys_SI.equivalent_dims(dim, units.mass / units.volume)

    # SymbolIndexedNew
    expr = clone_as_indexed(symbols.energy)
    dim = collect_dimension(expr)
    assert dimsys_SI.equivalent_dims(dim, units.energy)


def test_rest() -> None:
    expr = sympify(4)
    dim = collect_dimension(expr)
    assert dimsys_SI.is_dimensionless(dim)
