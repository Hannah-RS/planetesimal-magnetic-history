#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from parameters import rc, bpart

def gamma_aubert(f):
    """
    Scaling prefactor gamma (eqn 18) from p=gamma*Raq from Aubert (2009)

    Parameters
    ----------
    f : float
        fractional inner core radius.

    Returns
    -------
    gamma

    """
    # have simplifiedusing ri = f8rc and cancelled thing
    # have broken gamma into chunks to avoid typos
    a= (3*(1-f)**2)/(2*rc**2*(1-f**3))
    b = bpart*((3*rc**2*(1-f**5))/(5*(1-f**3))-f**2*rc**2)
    c = (1-bpart)*(rc**2-(3*rc**2*(1-f**5))/(5*(1-f**3)))
    
    gamma = a*(b+c)
    return gamma