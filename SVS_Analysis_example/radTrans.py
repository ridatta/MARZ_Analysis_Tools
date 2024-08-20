import numpy as np
from scipy.interpolate import interp1d
import pickle
from scipy.integrate import cumtrapz
from sklearn.neighbors import KNeighborsRegressor
from scipy import sparse
from scipy.sparse.linalg import spsolve







def simulate(ni0,sig,T0,E,smax,wl,fwhm=1.624):
    Iout = generateSpectrum(ni0,sig,T0,E,smax)
    Iout = convertToExp(wl,E,Iout,fwhm)
    Iout = (Iout - Iout.min()) / (Iout.max() - Iout.min()) # normalize between [0,1]
    return Iout

def simulateV2(ni0,sig,T0,delT,E,smax,wl,fwhm=1.624):
    Iout = generateSpectrumV2(ni0,sig,T0,delT,E,smax)
    Iout = convertToExp(wl,E,Iout)
#     Iout = (Iout - Iout.min()) / (Iout.max() - Iout.min()) # normalize between [0,1]
    return Iout

def simulateV3(N,sig,T,E,smax,wl,fwhm=1.624):
    Iout = np.zeros(E.shape)
    for ii in range(N.shape[0]):
        Iout = Iout + generateSpectrumV3(N[ii],sig,T[ii],E,smax)
    
    Iout = convertToExp(wl,E,Iout)
    Iout = (Iout - Iout.min()) / (Iout.max() - Iout.min()) # normalize between [0,1]

    return Iout

def baseline_als(y, lam, p, niter=10):
  L = len(y)
  D = sparse.diags([1,-2,1],[0,-1,-2], shape=(L,L-2))
  w = np.ones(L)
  for i in range(niter):
    W = sparse.spdiags(w, 0, L, L)
    Z = W + lam * D.dot(D.transpose())
    z = spsolve(Z, w*y)
    w = p * (y > z) + (1-p) * (y < z)
  return z

def gaussian(x,mu,sig): # Generates a gaussian distribution for a given mean mu and standard deviation sig
    y = np.exp(-(x-mu)**2 / (2 * sig**2))
    return y / y.max()

def add_padding(x):
    y1 = x[0] * np.ones(20)
    y2 = x[-1] * np.ones(20)
    arr = np.hstack([y1,x,y2])
    return arr

def depad(x):
    x = x[20:-20]
    return x

def convertToExp(wl,E,I,fwhm=1.624):
    # phototon energy from sim
    # Simulated spectrum
    
    # convert to nm space
    lam = 1239.8 / E # [nm]
    
    # Match dimensionality
    inp = interp1d(lam,I);
    Iout = inp(wl)
    
    # convolve
    sig_instrum = fwhm / 2.355
    s1 = gaussian(wl,wl[0]+(wl[-1]-wl[0])/2,sig_instrum)
    
    I_conv = np.convolve(add_padding(Iout),add_padding(s1),'same')
    I_conv = depad(I_conv)
    I_conv = I_conv / I_conv.max()
    
    return I_conv

def getEmiOpa(n,T,n_list,T_list,inDir):
    # returns photon energy [eV], emissivity [W/cm^3/sr/eV], and opacity [cm^-1]
    nidx = np.where(n_list == n)
    Tidx = np.where(T_list == T)
    fname = "spect_%04d_%04d.ppd" % (nidx[0]+1,Tidx[0]+1)

    data = np.loadtxt(inDir + fname,skiprows=18)
    E = data[:,0] # eV
    emi = data[:,1] # erg/s/sr/cm^3/eV
    opa = data[:,2] # cm^2/g
    
    # Interpolate so there is a fixed no. of spectral bins
    emi_inter = interp1d(E,emi); opa_inter = interp1d(E,opa)
    Enew = np.linspace(E.min(),E.max(),400)
    emi_new =  emi_inter(Enew); opa_new =  opa_inter(Enew)
    
    # Convert units
    A = 27.0
    rho = A * 1.6726e-24 * n # gm/cm^3
    opa_conv = opa_new * rho # cm^-1
    
    emi_conv = emi_new / 1e7 # W/cm^3/sr/eV
    
    return Enew,emi_conv,opa_conv

def gaussian(x,mu,sig):
    y = np.exp(-(x-mu)**2 / (2 * sig**2))
    return y / y.max()

def scale_ni(x):
    data = np.loadtxt("./Data/valuesV3.csv",
                 delimiter=",")
    n_list = data[:,0] # cm^-3
    T_list = data[:,1] # eV
    return (np.log10(x)-np.log10(n_list.min())) / (np.log10(n_list.max()) - np.log10(n_list.min()))

def scale_T(x):
    data = np.loadtxt("./Data/valuesV3.csv",
                 delimiter=",")
    n_list = data[:,0] # cm^-3
    T_list = data[:,1] # eV
    return (x-T_list.min()) / (T_list.max()-T_list.min())

def descale_ni(x):
    data = np.loadtxt("./Data/valuesV3.csv",
                 delimiter=",")
    n_list = data[:,0] # cm^-3
    T_list = data[:,1] # eV
    return np.around(10**((np.log10(n_list.max()) - np.log10(n_list.min())) * x + np.log10(n_list.min())),decimals=3)

def descale_T(x):
    data = np.loadtxt("./Data/valuesV3.csv",
                 delimiter=",")
    n_list = data[:,0] # cm^-3
    T_list = data[:,1] # eV
    return np.round(x * (T_list.max()-T_list.min()) + T_list.min(),decimals=2)

def getEmi(ni,Te): # Returns emissivity [W/cm^3/sr/eV] for given ion density [cm^-3] and temperature [eV]
    ni = scale_ni(ni)
    Te = scale_T(Te)
    knr_emi = pickle.load(open('./models/knr_emi.sav', 'rb'))  
    emi = knr_emi.predict(np.array([[ni,Te]])) # W/cm^3/sr/eV
    return emi[0]

def getOpa(ni,Te): # Returns opacity [cm^-1] for given ion density [cm^-3] and temperature [eV]
    ni = scale_ni(ni)
    Te = scale_T(Te)
    knr_opa = pickle.load(open('./models/knr_opa.sav', 'rb'))  
    opa = knr_opa.predict(np.array([[ni,Te]])) # cm^-1
    return opa[0]


def cumtrapz2(y,s,axis):
    Y = cumtrapz(y,x=s,axis=axis)
    Y = np.vstack([np.zeros(y.shape[1]),Y])
    return Y
    
def radtran(s,Iin,absr,emi): # Solves 1-D radaition transport for arbitraty emisisvity and opacity along s
    # s [m] = 1 x m vector between x = 0 and x = x_out
    # I_in [W/m^2/sr/eV] = 1 x n vector;input intensity spectrum descretized over n energy
    # datapoints
    # absr [m^-1] = m x n vector; rows correspond to spectral values; columns =
    # spatial values
    # emi [W/m^3/sr/eV] = m x n vector; rows correspond to spectral values; columns =
    # spatial dependendant values
    Iout = Iin * np.exp(-1 * np.trapz(absr,x=s,axis=0)) + np.trapz(emi * np.exp(-1 * cumtrapz2(absr,s,axis=0)),x=s,axis=0)
    return Iout

def radtranV2(s,Iin,absr,emi): # Solves 1-D radaition transport for arbitraty emisisvity and opacity along s
    # s [m] = 1 x m vector between x = 0 and x = x_out
    # I_in [W/m^2/sr/eV] = 1 x n vector;input intensity spectrum descretized over n energy
    # datapoints
    # absr [m^-1] = m x n vector; rows correspond to spectral values; columns =
    # spatial values
    # emi [W/m^3/sr/eV] = m x n vector; rows correspond to spectral values; columns =
    # spatial dependendant values
    tau_pr = cumtrapz2(absr,s,axis=0)
    tau =tau_pr[-1]
    S = emi / absr
    Iout = Iin * np.exp(-1 * tau) + np.trapz(S * np.exp(-1 * tau_pr),x=tau_pr,axis=0)
    return Iout

def simulate_combined(n,T,sig,E,wl,fwhm=1.624):
    # n, T, sig = vectors to capture vatiation across LOS

    out = np.zeros(E.shape)
    for ii in range(n.shape[0]):
        out = out + generateSpectrum(n[ii],sig[ii],T[ii],E,smax=2*2.355*sig[ii],normalize=False) # add the spectral intenisty, W/m^2/sr/eV
    out = convertToExp(wl,E,out,fwhm=1.624)
    return out / out.max() 


def generateSpectrum(ni0,sig,T0,E,smax=10e-3,normalize=True): 
# generates output intensity vs wavelength  for gaussian or uniform density distributions, uniform temp. after solving raditaion transport in 1-D
# ni0 = peak denisty [cm^-3]
# T0 = Temperature [eV]
# sig = std of gaussian densty distribution [m], set to 0 for unifomr density
# E = photon energy [eV]
# normalize: if set to TRUE, the outpus intensity is divided by max value

    # (1) density and temp
    m = 30 # no. of space steps
    s = np.linspace(0,smax,m) # [m], space
    
    if (sig == 0):
        ni = ni0 * np.ones(m)
    else:
        ni = ni0* gaussian(s,s[-1]/2,sig) # cm^-3
    Te = T0 * np.ones(m) # eV
    ni[np.where(ni<1e16)] = 1e16 # min. density
    
    # (2) Emisivity and Opacity
    n = E.shape[0] # No. spectral bins
    EM = np.zeros((m,n))
    OP = np.zeros((m,n))
    
    for ii in range(m):
        emi = getEmi(ni[ii],Te[ii])
        opa = getOpa(ni[ii],Te[ii]) 
        EM[ii,:] = emi
        OP[ii,:] = opa
    
    # (3) Radiation Transport
    Iin = np.zeros(n) # W/m^2/sr/eV
    EM_si = EM * 1e6 # W/m^3/sr/eV
    OP_si = OP * 1e2 # m^-1

    Iout = radtran(s,Iin,OP_si,EM_si) # W/m^2/sr/eV

    if normalize:
        Iout = Iout / Iout.max()

    return Iout


