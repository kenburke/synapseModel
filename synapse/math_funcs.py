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
    
def diffuse(r, i_ca, D_ca, k_on, k_d, B_total, Ca_basal):
    """
    Linearized Buffer Approximation
    Approximates Ca concentration as a function of distance from point source
    
    INPUT:
    - numpy array of distances in nanometers
    - instantaneous current through channel in picoamps
    - Various constants related to diffusion / buffering:
        - Diffusion Coefficient of free calcium (microns^2/s)
        - Rate of calcium binding to free buffer (M^-1 s^-1)
        - Dissociation constant (micromolar)
        - Total buffer concentration (molar)
        - basal calcium concentration (micromolar)
        
    OUTPUT:
    - numpy array of concentrations in micromolar (corresponding to input distances)
    
    """
    
    #Faraday's constant
    F = 96485.33289
    adjusted_ica = i_ca * 1e-12 # into amps
    adjusted_r = r / 1e3 # into microns

    B_free = (B_total * k_d) / (k_d + Ca_basal)
    
    lamb = (D_ca / (k_on * B_free)) ** 0.5
    conc = (adjusted_ica / (4.0 * np.pi * F * D_ca * adjusted_r)) * np.exp(-adjusted_r/lamb)
    
    conc *= 1e21
    
    return conc
    