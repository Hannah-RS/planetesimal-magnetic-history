#!/usr/bin/env python3
# -*- coding: utf-8 -*-
def Ephi(f):
    """
    Entropy associated with Ohmic dissipation caused by the dynamo. Based on N09
    

    Parameters
    ----------
    f : float
        fractional inner core radius

    Returns
    -------
    float,
     entropy due to Ohmic dissipation

    """
    from parameters import rc, phiv
    import numpy as np
    
    Voc = 4/3*np.pi*(1-f**3)*rc # volume of the outer core
    return Voc*phiv