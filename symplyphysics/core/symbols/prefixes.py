from collections import namedtuple

Prefixes = namedtuple("Prefixes", [
    "yotta", "zetta", "exa", "peta", "tera", "giga", "mega", "kilo", "hecto", "deca", "deci",
    "centi", "milli", "micro", "nano", "pico", "femto", "atto", "zepto", "yocto"
])

prefixes = Prefixes(
    yotta=10**24,
    zetta=10**21,
    exa=10**18,
    peta=10**15,
    tera=10**12,
    giga=10**9,
    mega=10**6,
    kilo=10**3,
    hecto=10**2,
    deca=10**1,
    deci=10**-1,
    centi=10**-2,
    milli=10**-3,
    micro=10**-6,
    nano=10**-9,
    pico=10**-12,
    femto=10**-15,
    atto=10**-18,
    zepto=10**-21,
    yocto=10**-24,
)
