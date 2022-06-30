#!/usr/bin/env python3
# -*- coding: utf-8 -*-
def stop(t,Tc,Ek=True,Er=True,Es=True,Eg=True,El=True):
    """
    Event to trigger when entire core has solidified to stop the integration

    Parameters
    ----------
    t : float
        time, s
    Tc : float
        core temperature

    Returns
    -------
    f-1 - zero crossing when f=1
    """
    from f_def import f_calc
    from parameters import Tc0, f0, B
    f=f_calc(Tc,Tc0,f0,B)
    return f-1