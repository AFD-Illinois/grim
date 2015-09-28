# Parallel movie rendering code:
# $ mpirun -np 8 python movie.py
#
# Code will spit out .png images. To make a movie do:
#
# ffmpeg -f image2 -i rho%04d.png -vcodec mpeg4 -mbd rd -trellis 2 -cmp 2 -g 300 -pass 1 -r 25 -b 18000000 movie.mp4
#
# [Replace ffmpeg by avconv wherever available]


import h5py
import glob
import numpy as np
import pylab as pl
from scipy.integrate import trapz
from mpl_toolkits.axes_grid1 import make_axes_locatable
import yt
yt.enable_parallelism()

# Set plot parameters to make beautiful plots
pl.rcParams['figure.figsize']  = 12, 7.5
pl.rcParams['lines.linewidth'] = 2.
pl.rcParams['font.family']     = 'serif'
pl.rcParams['font.weight']     = 'bold'
pl.rcParams['font.size']       = 30 
pl.rcParams['font.sans-serif'] = 'serif'
pl.rcParams['text.usetex']     = True
pl.rcParams['axes.linewidth']  = 1.5
pl.rcParams['axes.titlesize']  = 'medium'
pl.rcParams['axes.labelsize']  = 'medium'

pl.rcParams['xtick.major.size'] = 8
pl.rcParams['xtick.major.width'] = 2
pl.rcParams['xtick.minor.size'] = 4    
pl.rcParams['xtick.minor.width'] = 2
pl.rcParams['xtick.major.pad']  = 8    
pl.rcParams['xtick.minor.pad']  = 8    
pl.rcParams['xtick.color']      = 'k'    
pl.rcParams['xtick.labelsize']  = 'medium'
pl.rcParams['xtick.direction']  = 'in'   

pl.rcParams['ytick.major.size'] = 8    
pl.rcParams['ytick.major.width'] = 2
pl.rcParams['ytick.minor.size'] = 4    
pl.rcParams['ytick.minor.width'] = 2
pl.rcParams['ytick.major.pad']  = 8    
pl.rcParams['ytick.minor.pad']  = 8    
pl.rcParams['ytick.color']      = 'k'    
pl.rcParams['ytick.labelsize']  = 'medium'
pl.rcParams['ytick.direction']  = 'in'    

filepath = '/home/astro/work/frames_for_torus_movie/VaC-ModMirror_b/'
frames_folder_path = '/home/astro/work/frames_for_torus_movie/frames/'
file_prefix = 'plot_torus_'

# Load the primitive variables
data_files = np.sort(glob.glob(filepath+'/data*.h5'))[-1000:]


N1 = 280
N2 = 256
NG = 3
A_SPIN = 0.9375
H_SLOPE = 0.3
R_IN = 0.98*(1. + np.sqrt(1. - A_SPIN**2.))
R_OUT = 63.
X1_START = np.log(R_IN)
X2_START = 1e-20
DX1 = (np.log((R_OUT)/(R_IN))/N1)
DX2 = ((1. - 2.*X2_START)/N2)
i = np.arange(-NG, N1+NG, 1)
j = np.arange(-NG, N2+NG, 1)
X1 = (X1_START + (i + 0.5)*DX1)
X2 = (X2_START + (j + 0.5)*DX2)
R = np.exp(X1)
THETA = np.pi*(X2) + ((1 - H_SLOPE)/2.)*np.sin(2.*np.pi*(X2))

# Load the geometry
gCov = np.zeros([N2, N1, 4, 4])
gCon = np.zeros([N2, N1, 4, 4])
gammaDownDownDown = np.zeros([N2, N1, 4, 4, 4])
gammaUpDownDown = np.zeros([N2, N1, 4, 4, 4])

dumpfile = filepath + 'gcov.h5'
datafile = h5py.File(dumpfile, "r")
for mu in xrange(4):
    for nu in xrange(4):
        gCov[:, :, mu, nu] = datafile[datafile.items()[0][0]][:][:, :, nu + 4*mu]
        
dumpfile = filepath + 'gcon.h5'
datafile = h5py.File(dumpfile, "r")
for mu in xrange(4):
    for nu in xrange(4):
        gCon[:, :, mu, nu] = datafile[datafile.items()[0][0]][:][:, :, nu + 4*mu]
        
dumpfile = filepath + 'gammaUpdowndown.h5'
datafile = h5py.File(dumpfile, "r")
for eta in xrange(4):
    for mu in xrange(4):
        for nu in xrange(4):
            gammaUpDownDown[:, :, nu, mu, eta] = \
                datafile[datafile.items()[0][0]][:][:, :, nu + 4*(mu + 4*eta)]
                
g = np.zeros([N2, N1])
for j in xrange(N2):
    for i in xrange(N1):
        g[j, i] = np.sqrt(-np.linalg.det(gCov[j, i, :, :]))

def returnFluidElement(data):
    
    rho = data[:, :, 0]
    u   = data[:, :, 1]
    adiabatic_index = 4./3
    P   = (adiabatic_index - 1.)*u
    cs  = np.sqrt(adiabatic_index*P/(rho + adiabatic_index * u) )
    T   = P/rho

    uGrid = data[:, :, 2:5]
    bGrid = data[:, :, 5:8]

    gamma = 1. + np.zeros([N2, N1])
    alpha = np.sqrt(-gCon[:,:,0,0])**(-1)
    
    for i in xrange(0, 3):
        for j in xrange(0, 3):
            gamma += gCov[:,:,i+1,j+1]*uGrid[:,:,i]*uGrid[:,:,j]
            
    gamma = np.sqrt(gamma)
    uCon = np.zeros([N2, N1, 4])
    uCon[:,:,0] = gamma / alpha
    for i in xrange(1, 4):
        uCon[:,:,i] = uGrid[:,:,i-1]-gamma*alpha*gCon[:,:,0,i]
        
    uCov = np.zeros([N2, N1, 4])
    for mu in xrange(0, 4):
        for nu in xrange(0, 4):
            uCov[:,:,mu]+=gCov[:,:,mu,nu]*uCon[:,:,nu]
            
    bCon = np.zeros([N2, N1, 4])    
    for i in xrange(0, 3):
        bCon[:,:,0]+=bGrid[:,:,i]*uCov[:,:,i+1]
        
    for i in xrange(0, 3):
        bCon[:,:,i+1]=(bGrid[:,:,i]+bCon[:,:,0]*uCon[:,:,i+1])/uCon[:,:,0]
        
    bCov = np.zeros([N2, N1, 4])
    for mu in xrange(0, 4):
        for nu in xrange(0, 4):
            bCov[:,:,mu]+=gCov[:,:,mu,nu]*bCon[:,:,nu]
            
    bSqr = np.zeros([N2, N1])
    for mu in xrange(0, 4):
        bSqr += bCov[:,:,mu]*bCon[:,:,mu]
        
    ViscousCoeff = 1.
    ConductionCoeff = 1.
    betaV = 0.5/ViscousCoeff/cs**2/rho
    betaC = 1./ConductionCoeff/cs**2/rho/T
    deltaP = data[:, :, -1]*np.sqrt(T/betaV)
    q      = data[:, :, -2]*np.sqrt(T/betaC)

    B1 = data[:, :, 5]
    B2 = data[:, :, 6]
    B3 = data[:, :, 7]

    return {'rho': rho, 'u' : u , 'P' : P, 'cs' : cs , 'T' : T, \
            'uCon' : uCon, 'uCov' : uCov, 'bCon' : bCon, 'bCov' : bCov, \
            'bSqr' : bSqr, 'q' : q, 'deltaP' : deltaP, 'B1' : B1, 'B2' : B2, 'B3' : B3}

def returnMagneticVectorPotential(elem):

    B1 = elem['B1']
    B2 = elem['B2']
    vectorPotential = np.zeros([N2, N1])
    N_start = 1
    for j in xrange(N2):
        for k in xrange(N1):
           vectorPotential[j, k] = trapz(g[j, :k]*B2[j, :k], dx=DX1) - trapz(g[:j, k]*B1[:j, k], dx=DX2)
            
    return vectorPotential

XMin = 0
XMax = 25
YMin = -25
YMax = 25

R_grid, THETA_grid = np.meshgrid(R[NG:-NG], THETA[NG:-NG])
X = R_grid*np.sin(THETA_grid)
Y = R_grid*np.cos(THETA_grid)
alpha = np.sqrt(-gCon[:,:,0,0])**(-1)
PolGamma = 4./3

minorLocatorX   = pl.FixedLocator(np.linspace(0,60,13))
minorLocatorY   = pl.FixedLocator(np.linspace(-20,20,9))
majorLocatorX   = pl.FixedLocator(np.linspace(0,60,7))
majorLocatorY   = pl.FixedLocator(np.linspace(-20,20,5))

for file_number, dump_file in yt.parallel_objects(enumerate(data_files)):

    print "File number = ", file_number
    frame_index = file_number

    data_file = h5py.File(dump_file, "r")
    primVars  = data_file['primVars']

    elem = returnFluidElement(primVars)
    A_plot = returnMagneticVectorPotential(elem)
    N_start = 1
            
    fig, axes = pl.subplots(nrows=1, ncols=3)
        
    fig.set_size_inches((24, 12))
    
    fontsize_ticks = 20
    
    pl.subplot(131)
    ax = pl.gca()
    ax.patch.set_facecolor('black')

    ax.xaxis.set_major_locator(majorLocatorX)
    ax.xaxis.set_minor_locator(minorLocatorX)
    ax.xaxis.set_ticklabels([0,10,20,30,40,50],fontsize=fontsize_ticks)
    ax.yaxis.set_major_locator(majorLocatorY)
    ax.yaxis.set_minor_locator(minorLocatorY)
    ax.yaxis.set_ticklabels([-20,-10,0,10,20],fontsize=fontsize_ticks)
    pl.title('$\\log_{10}(\\rho)$')

    pl.axis([XMin, XMax, YMin, YMax])
    ax.set_aspect('equal')

    colouraxis = np.linspace(-6., 0., 100)
    im=pl.contourf(X, Y, pl.log10(elem['rho']), colouraxis, cmap='YlOrRd')
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.1)
    cbar = fig.colorbar(im,ticks=[-6,-5,-4,-3,-2,-1,0], cax=cax)
    
    pl.subplot(132)
    ax = pl.gca()
    ax.patch.set_facecolor('black')

    ax.xaxis.set_major_locator(majorLocatorX)
    ax.xaxis.set_minor_locator(minorLocatorX)
    ax.xaxis.set_ticklabels([0,10,20,30,40,50],fontsize=fontsize_ticks)
    ax.yaxis.set_major_locator(majorLocatorY)
    ax.yaxis.set_minor_locator(minorLocatorY)
    ax.yaxis.set_ticklabels([-20,-10,0,10,20],fontsize=fontsize_ticks)
    pl.title('$\\log_{10}(q/\\rho c_s^3)$')
    
    pl.axis([XMin, XMax, YMin, YMax])
    ax.set_aspect('equal')
    
    colouraxis = np.linspace(-3, 0., 100)
    im=pl.contourf(X, Y, np.log10(abs(elem['q']/elem['rho']/elem['cs']**3.)), \
                   colouraxis, extend='min', cmap='gist_stern_r')

    vectorPotentialContourLines = np.linspace(-.22, .13, 100)
    pl.contour(X[N_start:-N_start, N_start:-N_start],
        Y[N_start:-N_start, N_start:-N_start],
        A_plot[N_start:-N_start, N_start:-N_start],
        vectorPotentialContourLines, colors='black', linestyles='solid',
        alpha=0.25)

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.1)
    cbar = fig.colorbar(im,ticks=[-4, -3,-2,-1,0], cax=cax)
        
    
    pl.subplot(133)
    ax = pl.gca()
    ax.patch.set_facecolor('black')

    ax.xaxis.set_major_locator(majorLocatorX)
    ax.xaxis.set_minor_locator(minorLocatorX)
    ax.xaxis.set_ticklabels([0,10,20,30,40,50],fontsize=fontsize_ticks)
    ax.yaxis.set_major_locator(majorLocatorY)
    ax.yaxis.set_minor_locator(minorLocatorY)
    ax.yaxis.set_ticklabels([-20,-10,0,10,20],fontsize=fontsize_ticks)
    pl.title('$\\Delta P/b^2$')
    
    pl.axis([XMin, XMax, YMin, YMax])
    ax.set_aspect('equal')

    colouraxis = np.linspace(-1.1, 0.6, 100)
    im=pl.contourf(X, Y, elem['deltaP']/elem['bSqr'], colouraxis, cmap='bwr')
    
    pl.contour(X[N_start:-N_start, N_start:-N_start],
        Y[N_start:-N_start, N_start:-N_start],
        A_plot[N_start:-N_start, N_start:-N_start],
        vectorPotentialContourLines, colors='black', linestyles='solid',
        alpha=0.25)
    
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.1)
    cbar = fig.colorbar(im,ticks=[-1,-0.5,0,0.5], cax=cax)
    
    pl.tight_layout()
    pl.savefig(frames_folder_path + '/' + file_prefix + '%04d'%frame_index + \
            '.png')
    pl.close()

    data_file.close()
