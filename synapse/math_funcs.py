import numpy as np

def hill(ca,S=1,ec50=0.7,n=3.72):
    """

    Hill function

    INPUT:
    - numpy array of [Ca]
    - Hill parameters (S, EC_50, N)

    OUTPUT:
    - numpy array of size ca that is Hill(ca)

    """
    pr = (S * (ca ** n)) / (ec50 ** n + ca ** n)
    return pr
    
def ampa(time,quantal_size=-5,tau=(1/1000,5/1000)):
    """
    
    AMPA kinetics as difference of exponentials
    
    INPUT:
    - numpy array of time
    - quantal size (maximal amplitude, pA)
    - time constants tau (in ms, tuple of (fast,slow))
    
    OUTPUT:
    - numpy array of size time
    
    """
    i_ampa = np.exp(-time/tau[1]) - np.exp(-time/tau[0])
    i_ampa *= quantal_size/max(i_ampa)
    
    return i_ampa