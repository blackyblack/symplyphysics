def filter_zeroes(list_):
    return list(filter(lambda e: e != 0, list_))


# list_ contains map as an element
def filter_map_zeroes(key_, list_):
    return list(filter(lambda e: e[key_] != 0, list_))


def filter_negative(list_):
    return list(filter(lambda e: e >= 0, list_))


# list_ contains map as an element
def filter_map_negative(key_, list_):
    return list(filter(lambda e: e[key_] >= 0, list_))
