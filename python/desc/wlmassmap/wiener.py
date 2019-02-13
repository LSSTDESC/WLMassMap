from vmad import Builder, autooperator
from vmad.lib import fastpm, mpi, linalg
import wiener_cosmo2d as cosmo2d
from abopt.algs.trustregion import TrustRegionCG
import numpy 
from abopt.algs.lbfgs import LBFGS
from abopt.algs.lbfgs import pre_scaled_direct_bfgs, pre_scaled_inverse_dfp
from abopt.linesearch import minpack, backtrace, exact
from abopt.abopt2 import real_vector_space, Problem, VectorSpace
from nbodykit.cosmology import Planck15, LinearPower
import warnings



oldprint = print

noise_file = "./LBFGS-bar-truth.inoise_variance.npy"
data_file = "./LBFGS-bar-truth.d.npy"


inoise_variance = numpy.load(noise_file)
data =  numpy.load(data_file)

def flat_WF_map(gmap,nmap, BoxSize=10. ,smooth_scale=0.1,maxiter = 1000,fid_ps = Planck15,verbose = False):
    """
    Wiener Filter Map Making: As in Horowitz, Seljak, Aslanyan (2018)

    nmap : ndarray
       Inverse variance map of the data. A noisemap with same shape as gamma map. 
    BoxSize : float
       Sidelength of box in Mpc; only square boxes work (it is likely easy to fix or just pad input)
    smooth_scale: float
       Smoothing scale on reconstructed map to maintain numerical accuracy
    maxiter: int
       Maximum number of iterations for reconstruction
    fid_ps: function
       Fiducial PS to start optimization around. This shouldn't matter except to speed up 
       optimization (check manually WMAP9 vs. Planck15 if curious)...
    verbose: bolean
       additional information is outputted, including numpy arrays at each iteration step

    """

    ## Checking inputs... Errors
    if gmap.shape[0] != gmap.shape[1]:
        raise NotImplementedError("Only allow square input shear maps! Got map of size...", gmap.shape[0],gmap.shape[1])
    if gmap.shape[2] != 2:
        raise NotImplementedError("Require two shear maps, got...",gmap.shape[2])
    if gmap.shape != nmap.shape:
        raise NotImplementedError("Shear map must have same shape as noise map!")
    if gmap.shape[0] < 10:
        raise NotImplementedError("Shear map has too small pixel size;", gmap.shape[0])
    
    ## Checking inputs... Warnings
    if smooth_scale >0.5:
        warnings.warn("Smoothing scale is very large; "+str(smooth_scale)+"! View results cautiously!" , UserWarning)

    checksum = abs(data[:,0,0])+abs(data[:,0,1])+abs(data[0,:,0])+abs(data[0,:,1])
    checksum += abs(data[:,-1,0])+abs(data[:,-1,1])+abs(data[-1,:,0])+abs(data[-1,:,1])
    if numpy.sum(checksum) > numpy.abs(data).mean():
        warnings.warn("This Wiener Filter implementation assumes periodic boundary conditions. There are non-zero boundary elements in the input shear maps; view results cautiously! Pad input if need be!" , UserWarning)

    ## Calling main function
    bf_signal,bf_gamma = __flat_WF_map(gmap,nmap,BoxSize=BoxSize,smooth_scale=smooth_scale,maxiter=maxiter,fid_ps=fid_ps,verbose=verbose)

    return bf_signal

def __flat_WF_map(gmap,nmap,BoxSize=10. ,smooth_scale=0.1,maxiter = 1000,fid_ps = Planck15,verbose = False):


    Nsize = data.shape[0]

    pm = cosmo2d.ParticleMesh([Nsize, Nsize], BoxSize=BoxSize)

    try:
        powerspectrum = LinearPower(fid_ps, 0)
    except:
        powerspectrum = fid_ps

    wn = pm.generate_whitenoise(556, unitary=True)

    x = wn[...]
    x = numpy.stack([x.real, x.imag], -1)

    ForwardModelHyperParameters = dict(
                powerspectrum=powerspectrum,
                pm=pm,
                smooth_scale=smooth_scale)

    ForwardOperator = cosmo2d.FastPMOperator.bind(**ForwardModelHyperParameters)


    problem = cosmo2d.ChiSquareProblem(pm.comm,
            ForwardOperator,
            [
                cosmo2d.LNResidualOperator,
                cosmo2d.NLResidualOperator.bind(d=data, isigma=inoise_variance),
            ]
            )


    problem.maxradius = 100
    problem.initradus = 1
    problem.cg_rtol = 0.1
    problem.cg_maxiter= 10

    def monitor(state):
        problem.save('./output/lb-bar-%04d' % state['nit'], state)
        print(state)
        #print("Time: ",time.time()-start_time)

    lbfgs = LBFGS(maxiter=maxiter, linesearch=exact, diag_update=pre_scaled_direct_bfgs)
    if verbose:
        print('objective(truth) =', problem.f(x), 'expecting', pm.Nmesh.prod() * len(problem.residuals))
        x1 = lbfgs.minimize(problem, x * 0.001, monitor=monitor)
    else:
        x1 = lbfgs.minimize(problem, x * 0.001)

    s_sol,fs_sol = problem.state(x1)

    return np.fft.irfft2(s_sol*Nmesh*Nmesh),fs_sol
