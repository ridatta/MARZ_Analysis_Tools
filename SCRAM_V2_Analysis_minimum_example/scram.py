import numpy as np
import pickle
from scipy.integrate import cumtrapz

def simulate_gauss(E,ni_h,T_h,ni_l,T_l,L,d,pos,m=800):
    n = E.shape[0]
    s = np.linspace(0,L,m) # m
    
    ni = ni_l * np.ones(s.shape) # cm^-3
    T = T_l * np.ones(s.shape) # eV
    
    
    for ii in range(pos.shape[0]):
        ni = ni +  (ni_h-ni_l)*gaussian(s,pos[ii],d/2.355) 
        T = T +  (T_h-T_l)*gaussian(s,pos[ii],d/2.355) 
    
    
    # Show emissivity and opacity vs photon energy vs s
    EM = np.zeros((m,n))
    OP = np.zeros((m,n))

    for ii in range(m):
        emi = getEmi(ni[ii],T[ii]) * 1e6 # W/m^3/sr/eV
        opa = getOpa(ni[ii],T[ii]) * 100 # m^-1
        EM[ii,:] = emi 
        OP[ii,:] = opa
        
    Iin = np.zeros(E.shape) # W/m^2/sr/eV
    Iout = radtran(s,Iin,OP,EM)
    return Iout, s, ni, T

def simulate_rect(E,ni_h,T_h,ni_l,T_l,L,d,m=400):
    n = E.shape[0]
    s = np.linspace(0,L,m) # m
    
    ni = ni_l * np.ones(s.shape) # cm^-3
    T = T_l * np.ones(s.shape) # eV
    
    ni[np.where((s>=0.5*L-d/2) & (s<=0.5*L+d/2))] = ni_h
    T[np.where((s>=0.5*L-d/2) & (s<=0.5*L+d/2))] = T_h
    
    # Show emissivity and opacity vs photon energy vs s
    EM = np.zeros((m,n))
    OP = np.zeros((m,n))

    for ii in range(m):
        emi = getEmi(ni[ii],T[ii]) * 1e6 # W/m^3/sr/eV
        opa = getOpa(ni[ii],T[ii]) * 100 # m^-1
        EM[ii,:] = emi 
        OP[ii,:] = opa
        
    Iin = np.zeros(E.shape) # W/m^2/sr/eV
    Iout = radtran(s,Iin,OP,EM)
    return Iout, s, ni, T


def simulate_homogenous(E,ni,T,L,m=400):
    n = E.shape[0]
    
    s = np.linspace(0,L,m) # m
    ni = ni * np.ones(s.shape) # cm^-3
    T = T * np.ones(s.shape) # eV
    
    # Show emissivity and opacity vs photon energy vs s
    EM = np.zeros((m,n))
    OP = np.zeros((m,n))

    for ii in range(m):
        emi = getEmi(ni[ii],T[ii]) * 1e6 # W/m^3/sr/eV
        opa = getOpa(ni[ii],T[ii]) * 100 # m^-1
        EM[ii,:] = emi 
        OP[ii,:] = opa
        
    Iin = np.zeros(E.shape) # W/m^2/sr/eV
    Iout = radtran(s,Iin,OP,EM)
    return Iout
    

def instrument_response(E,fwhm):
    # returns a gaussin function with FWHM 
    sig = fwhm / 2.355
    return gaussian(E,E.min() + 0.5 * (E.max()-E.min()),sig)

def cumtrapz2(y,s,axis):
    Y = cumtrapz(y,x=s,axis=axis)
    Y = np.vstack([np.zeros(y.shape[1]),Y])
    return Y
    
def radtran(s,Iin,absr,emi):
    # s [m] = 1 x m vector between x = 0 and x = x_out
    # I_in [W/m^2/sr/eV] = 1 x n vector;input intensity spectrum descretized over n
    # datapoints
    # absr [m^-1] = m x n vector; rows correspond to spectral values; columns =
    # spatial values
    # emi [W/m^3/sr/eV] = m x n vector; rows correspond to spectral values; columns =
    # spatial dependendant values
    Iout = Iin * np.exp(-1 * np.trapz(absr,x=s,axis=0)) + np.trapz(emi * np.exp(-1 * cumtrapz2(absr,s,axis=0)),x=s,axis=0)
    return Iout

def scale_ni(x):
    n_list = np.loadtxt('./models/ni.txt') # cm^-3
    return (np.log10(x)-np.log10(n_list.min())) / (np.log10(n_list.max()) - np.log10(n_list.min()))

def scale_T(x):
    T_list = np.loadtxt('./models/Te.txt')
#     return (x-T_list.min()) / (T_list.max()-T_list.min()) # linear
    return (np.log10(x)-np.log10(T_list.min())) / (np.log10(T_list.max()) - np.log10(T_list.min()))

def descale_ni(x):
    n_list = np.loadtxt('./models/ni.txt') # cm^-3
    return np.around(10**((np.log10(n_list.max()) - np.log10(n_list.min())) * x + np.log10(n_list.min())),decimals=3)

def descale_T(x):
    T_list = np.loadtxt('./models/Te.txt')
#     return np.round(x * (T_list.max()-T_list.min()) + T_list.min(),decimals=2) # linear
    return np.around(10**((np.log10(T_list.max()) - np.log10(T_list.min())) * x + np.log10(T_list.min())),decimals=2)

def getEmi(ni,Te):
    ni = scale_ni(ni)
    Te = scale_T(Te)
    knr_emi = pickle.load(open('./models/knr_scram_emi.sav', 'rb'))  
    emi = knr_emi.predict(np.array([[ni,Te]])) # W/cm^3/sr/eV
    return emi[0]

def getOpa(ni,Te):
    ni = scale_ni(ni)
    Te = scale_T(Te)
    knr_opa = pickle.load(open('./models/knr_scram_opa.sav', 'rb'))  
    opa = knr_opa.predict(np.array([[ni,Te]])) # cm^-1
    return opa[0]

def gaussian(x,mu,sig):
    y = np.exp(-(x-mu)**2 / (2 * sig**2))
    return y / y.max()