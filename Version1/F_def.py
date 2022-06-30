#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Function F from Nimmo 2009
:param f: float, r_i=fr_c

:returns F: float
"""
def F_calc(f):
    return (1/5+2/15*f**5-(f**2)/3)*f/(1-f**3)