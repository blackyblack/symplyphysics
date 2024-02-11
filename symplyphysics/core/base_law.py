from typing import Any

from sympy import solve, Eq, S
from symplyphysics import Symbol, Quantity, convert_to, print_expression
from symplyphysics.core.quantity_decorator import _assert_expected_unit
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
        self.symbols = law.atoms(Symbol)

        variables_keys = list(str(symbol) for symbol in law.atoms(Symbol))
        for key, symbol in zip(variables_keys, self.symbols):
            self.__dict__[key[:-1]] = symbol    # key[:-1] for deleting index in end of string

    def __getattr__(self, name: str) -> Any:
        return self.__dict__[name]

    def __setattr__(self, name: str, value: Any) -> None:
        self.__dict__[name] = value

    def print_law(self) -> str:
        return print_expression(self.law)

    def calculate_symbol_value(
        self,
        variable: Symbol,
        dict_for_subs: dict[Symbol, Quantity],
        convert_to_type: type[float | int] | None = None
    ) -> Quantity | Probability | Fraction | float | int:
        self.__check_symbol_in_equation(variable)
        for symbol in dict_for_subs.keys():
            self.__check_symbol_in_equation(symbol)

        self.__validate_input(dict_for_subs)
        solved = solve(self.law, variable, dict=True)[0][variable]
        result_expr = solved.subs(dict_for_subs)
        result = Quantity(result_expr)
        self.__validate_output(variable, result)

        if convert_to_type:
            return convert_to_type(convert_to(result, S.One).evalf())
        return result

    def __check_symbol_in_equation(self, symbol: Symbol) -> None:
        if symbol not in self.symbols:
            raise ValueError(f"The specified variable {symbol} is not in the law.")

    def __validate_input(self, input_dict: dict[Symbol, Quantity]) -> None:
        for symbol, quantity in input_dict.items():
            _assert_expected_unit(quantity, symbol, symbol.name[:-1], self.calculate_symbol_value.__name__)

    def __validate_output(self, variable: Symbol, result: Quantity) -> None:
        _assert_expected_unit(result, variable, "return", self.calculate_symbol_value.__name__)
