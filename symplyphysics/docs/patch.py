"""
This module provides a function `patch_sympy_evaluate` for disabling `sympy` evaluation in certain
situations for documentation purposes.
"""

import ast

_IMPORT_NODE = ast.ImportFrom("symplyphysics.core.processors",
    names=[ast.alias(name="disable_sympy_evaluation"),
    ast.alias(name="reset_sympy_evaluation")],
    level=0)

_DISABLE_NODE = ast.Expr(ast.Call(ast.Name("disable_sympy_evaluation", ctx=ast.Load()), [], []))

_ENABLE_NODE = ast.Expr(ast.Call(ast.Name("reset_sympy_evaluation", ctx=ast.Load()), [], []))


def patch_sympy_evaluate(module: ast.Module) -> ast.Module:
    """
    Patches ``module`` so that `sympy` does not evaluate the expressions within it for
    documentation purposes, except for nodes where explicit evaluation has been documented.
    """

    module.body.insert(1, _IMPORT_NODE)

    current_member_idx = -1
    disabled_node_idxs: list[int] = []
    last_documented_node = 0

    for idx, stmt in enumerate(module.body):
        if isinstance(stmt, ast.FunctionDef):
            doc = ast.get_docstring(stmt)
            if doc is None:
                continue
            current_member_idx = idx
            last_documented_node = idx
            continue

        if isinstance(stmt, ast.Assign):
            for target in stmt.targets:
                try:
                    name = getattr(target, "id")
                    if str(name).startswith("_"):
                        continue
                except AttributeError:
                    continue
                current_member_idx = idx
            continue

        if (isinstance(stmt, ast.Expr) and current_member_idx >= 0 and
                isinstance(stmt.value, ast.Constant)):
            last_documented_node = idx

            s = str(stmt.value.value)
            # Skip members that either require evaluation or don't require autogeneration
            if (":laws:sympy-eval::" in s or
                (":laws:symbol::" not in s and ":laws:latex::" not in s)):
                continue

            disabled_node_idxs.append(current_member_idx)

    # Delete code unrelated to documentation
    module.body = module.body[0:last_documented_node + 1]

    offset = 0
    for node_idx in disabled_node_idxs:
        module.body.insert(node_idx + offset, _DISABLE_NODE)
        offset += 1
        # Enable eval after docstring
        module.body.insert(node_idx + offset + 1, _ENABLE_NODE)
        offset += 1

    ast.fix_missing_locations(module)
    return module


__all__ = ["patch_sympy_evaluate"]
