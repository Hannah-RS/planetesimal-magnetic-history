#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Fe-FeS liquidus - linear liquidus approximation for eutectic solidification at 32% S at 1234K
liquid iron melting at 1810 K 
see Kathryn's thesis pg 180
"""
def fe_fes_liquidus(Xs):
    """

    Parameters
    ----------
    Xs : float
        wt % sulfur

    Returns
    -------
    None.

    """
    return 1810-18*Xs