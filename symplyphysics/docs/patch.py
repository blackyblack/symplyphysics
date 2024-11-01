import ast


def patch_sympy_evaluate(content: ast.Module) -> ast.Module:
    import_node = ast.ImportFrom("symplyphysics.core.processors",
        names=[
        ast.alias(name="disable_sympy_evaluation"),
        ast.alias(name="reset_sympy_evaluation")
        ],
        level=0)
    disable_node = ast.Expr(ast.Call(ast.Name("disable_sympy_evaluation", ctx=ast.Load()), [], []))
    enable_node = ast.Expr(ast.Call(ast.Name("reset_sympy_evaluation", ctx=ast.Load()), [], []))

    content.body.insert(1, import_node)

    current_member_idx: int = -1
    disabled_nodes: list[int] = []
    last_documented_node: int = 0
    for idx, e in enumerate(content.body):
        if isinstance(e, ast.FunctionDef):
            doc = ast.get_docstring(e)
            if doc is None:
                continue
            current_member_idx = idx
            continue
        if isinstance(e, ast.Assign):
            for t in e.targets:
                try:
                    name = getattr(t, "id")
                    if str(name).startswith("_"):
                        continue
                except AttributeError:
                    continue
                current_member_idx = idx
            continue
        if isinstance(e, ast.Expr) and current_member_idx >= 0 and isinstance(
                e.value, ast.Constant):
            last_documented_node = idx
            if str(e.value.value).find(":laws:sympy-eval::") >= 0:
                continue
            # Skip members that don't require autogeneration
            if str(e.value.value).find(":laws:symbol::") < 0 and str(
                    e.value.value).find(":laws:latex::") < 0:
                continue
            disabled_nodes.append(current_member_idx)

    # Delete unrelated to documentation code
    content.body = content.body[0:last_documented_node + 1]

    offset = 0
    for n in disabled_nodes:
        content.body.insert(n + offset, disable_node)
        offset = offset + 1
        # enable eval after docstring
        content.body.insert(n + offset + 2, enable_node)
        offset = offset + 1

    ast.fix_missing_locations(content)
    return content
