from sympy import solve, Eq, S
from sympy.physics.units import Dimension

from symplyphysics import Symbol, Quantity, errors, dimensionless, convert_to, \
    print_expression
from symplyphysics.core.symbols.fraction import Fraction
from symplyphysics.core.symbols.probability import Probability


class BaseLaw:
    """
    Example:
    --------

    from symplyphysics.core.base_law import BaseLaw
    from symplyphysics.laws.thermodynamics import laplas_pressure as laplas_law

    law = BaseLaw(laplas_law.law)
    result = law.calculate_symbol_value(
        laplas_law.laplas_pressure,
        dict_quantities={
            laplas_law.radius_of_curvature: Quantity(0.05 * units.meters),
            laplas_law.surface_tension_of_the_liquid: Quantity(7 * units.newtons / units.meter)
        }
    )
    print(result, result.dimension)

    or
    --

    from symplyphysics.core.base_law import BaseLaw
    from symplyphysics.laws.thermodynamics import laplas_pressure as laplas_law

    law = BaseLaw(laplas_law.law)
    result = law.calculate_symbol_value(
        law.laplas_pressure,
        dict_quantities={
            law.radius_of_curvature: Quantity(0.05 * units.meters),
            law.surface_tension_of_the_liquid: Quantity(7 * units.newtons / units.meter)
        }
    )
    print(result, result.dimension)
    """

    def __init__(self, law: Eq) -> None:
        self.law = law
        self.symbols = self.law.atoms(Symbol)

        # NOTE: mypy raised "BaseLaw has no attribute..." if use this attributes
        variables_keys = list(str(symbol) for symbol in self.symbols)
        for key, symbol in zip(variables_keys, self.symbols):
            self.__dict__[key[:-1]] = symbol    # key[:-1] for deleting index in end of string

    def print_law(self) -> str:
        return print_expression(self.law)

    def calculate_symbol_value(
        self,
        variable: Symbol,
        dict_quantities: dict[Symbol, Quantity],
        convert_to_type: type[float | int] | None = None
    ) -> Quantity | Probability | Fraction | float | int:
        self.__check_symbol_in_equation(variable)
        for symbol in dict_quantities.keys():
            self.__check_symbol_in_equation(symbol)

        dict_validate = self.__get_validate_input_dict(dict_quantities)
        self.__validate_input(dict_validate)
        solved = solve(self.law, variable, dict=True)[0][variable]
        result_expr = solved.subs(dict_quantities)
        result = Quantity(result_expr)
        self.__validate_output(variable, result)

        if convert_to_type:
            return convert_to_type(convert_to(result, S.One).evalf())
        return result

    def __check_symbol_in_equation(self, symbol: Symbol) -> None:
        if symbol not in self.symbols:
            raise ValueError(f"The specified variable {symbol} is not in the law.")

    def __get_validate_input_dict(self, dict_of_quantities: dict[Symbol, Quantity]) -> dict[Dimension, Dimension]:
        quantity_type = lambda quantity: quantity.dimension if isinstance(quantity, Quantity) else dimensionless
        dict_validate = {symbol.dimension: quantity_type(quantity) for symbol, quantity in dict_of_quantities.items()}
        return dict_validate

    def __validate_input(self, validate_dict: dict[Dimension, Dimension]) -> None:
        for symbol_dimension, quantity_dimension in validate_dict.items():
            if isinstance(quantity_dimension, Dimension) and symbol_dimension != quantity_dimension:
                raise errors.UnitsError(f"Expected dimension: {symbol_dimension}, received: {quantity_dimension}.")

    def __validate_output(self, variable: Symbol, result: Quantity) -> None:
        result_type = lambda quantity: quantity.dimension if isinstance(quantity, Quantity) else dimensionless
        variable_dimension = variable.dimension
        result_dimension = result_type(result)
        if variable_dimension != result_type(result):
            raise errors.UnitsError(f"Expected dimension: {variable_dimension}, received: {result_dimension}.")
