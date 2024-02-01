import functools


# This decorator is used to mark tests that show unsupported usage of some class, usually fields.
# For example, using effects in field function is forbidden, so we mark it with @unsupported_usage.
# You're not supposed to use class this way, but we cannot enforce it by Python code.
# Use tests marked with `unsupported_usage` as a reference on how field behaves when passed incorrect input.
def unsupported_usage(func):

    @functools.wraps(func)
    def wrapper_(*args, **kwargs):
        return func(*args, **kwargs)

    return wrapper_
