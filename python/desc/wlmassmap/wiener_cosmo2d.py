from vmad import Builder, autooperator
from vmad.lib import fastpm, mpi, linalg

import numpy

from pmesh.pm import ParticleMesh

@autooperator
class FastPMOperator:
    ain = [('x', '*')]
    aout = [('wn', '*'), ('s', '*'), ('fs', '*')]

    def main(self, x, powerspectrum, pm,smooth_scale):
        wnk = fastpm.as_complex_field(x, pm)

        wn = fastpm.c2r(wnk)
        rholnk = fastpm.induce_correlation(wnk, powerspectrum, pm)



        rholnk = fastpm.apply_transfer(rholnk, smoothing(smooth_scale))
      #  q = pm.generate_uniform_particle_grid()
      #  layout = fastpm.decompose(q, pm)

        #phi = fastpm.apply_transfer(rholnk, fastpm.fourier_space_laplace)
#       rho = fastpm.c2r(rholnk)
#        rho = linalg.add(rho, 1.0)
        dx_c = fastpm.apply_transfer(rholnk, fourier_space_neg_gradient_lownoise(0,pm))
        dx2_c = fastpm.apply_transfer(dx_c, fourier_space_neg_gradient_lownoise(0,pm))

        dy_c = fastpm.apply_transfer(rholnk, fourier_space_neg_gradient_lownoise(1,pm))
        dy2_c = fastpm.apply_transfer(dy_c, fourier_space_neg_gradient_lownoise(1,pm))

        dxdy_c = fastpm.apply_transfer(dy_c, fourier_space_neg_gradient_lownoise(0,pm))

        g1_c = linalg.mul(dxdy_c,-2.0)
        g2_int = linalg.mul(dy2_c,-1.0)
        g2_c = linalg.add(dx2_c,g2_int)

        g1_r = fastpm.c2r(g1_c)#g1_c)
        g2_r =fastpm.c2r(g2_c)#g2_c)

       # g1 = fastpm.readout(g1_r, q, layout)
       # g2 = fastpm.readout(g2_r, q, layout)

        inter = [g1_r,g2_r]
        fs = linalg.stack(inter,axis=-1)

        return dict(fs=fs, s=rholnk, wn=wn)


def smoothing(scale):
    def tf(k):
        k2 = sum(ki ** 2 for ki in k)
        return numpy.exp(- 0.5 * k2 * scale ** 2)
    return tf

def fourier_space_neg_gradient_lownoise(dir,pm):
    def kernel(k):
            cellsize = (pm.BoxSize[dir] / pm.Nmesh[dir])
            w = k[dir] * cellsize

            a = 1 / (6.0 * cellsize) * (8 * numpy.sin(w) - numpy.sin(2 * w))
            # a is already zero at the nyquist to ensure field is real
            return -1.0 * (1j * a)
    return kernel

@autooperator
class NLResidualOperator:
    ain = [('wn', '*'), ('s', '*'), ('fs', '*')]
    aout = [('y', '*')]
    def main(self, wn, s, fs, d, isigma):
        r = linalg.add(fs, d * -1)
        r = linalg.mul(r, isigma)

        return dict(y = r)

@autooperator
class SmoothedNLResidualOperator:
    ain = [('wn', '*'), ('s', '*'), ('fs', '*')]
    aout = [('y', '*')]
    def main(self, wn, s, fs, d, isigma, scale):
        r = linalg.add(fs, d * -1)
        def tf(k):
            k2 = sum(ki ** 2 for ki in k)
            return numpy.exp(- 0.5 * k2 * scale ** 2)
        c = fastpm.r2c(r)
        c = fastpm.apply_transfer(c, tf)
        r = fastpm.c2r(c)
        r = linalg.mul(r, isigma)
        return dict(y = r)

@autooperator
class LNResidualOperator:
    ain = [('wn', '*'), ('s', '*'), ('fs', '*')]
    aout = [('y', '*')]
    def main(self, wn, s, fs):
        r = linalg.add(wn, 0)
        #fac = linalg.pow(wn.Nmesh.prod(), -0.5)
        fac = linalg.pow(1024*1024, -0.5)

        r = linalg.mul(r, fac)
        return dict(y = r)

@autooperator
class ChiSquareOperator:
    ain = [('x', '*')]
    aout = [('y', '*')]

    def main(self, x, comm):
        chi2 = linalg.sum(linalg.mul(x, x))
        chi2 = mpi.allreduce(chi2, comm)
        return dict(y = chi2)

from abopt.abopt2 import Problem as BaseProblem, VectorSpace

class ChiSquareProblem(BaseProblem):
    """ Defines a chisquare problem, which is

        .. math::

            y = \sum_i [R_i(F(s))]^2

        F : forward_operator

        [R_i] : residuals

        comm : The problem is defined on a MPI communicator -- must be consistent
        with forward_operator and residuals. (usually just pm.comm or MPI.COMM.WORLD)

    """
    def __init__(self, comm, forward_operator, residuals):
        self.residuals = residuals
        self.forward_operator = forward_operator
        self.comm = comm

        with Builder() as m:
            x = m.input('x')
            y = 0
            wn, s, fs = forward_operator(x)
            # fixme: need a way to directly include a subgraphs
            # rather than building it again.
            for operator in self.residuals:
                r = operator(wn, s, fs)
                chi2 = ChiSquareOperator(r, comm)
                y = linalg.add(y, chi2)
            m.output(y=y)

        def objective(x):
            return m.compute(vout='y', init=dict(x=x))

        def gradient(x):
            y, [vjp] = m.compute_with_vjp(init=dict(x=x), v=dict(_y=1.0))
            return vjp

        def hessian_vector_product(x, v):
            Dv = 0

            replay = forward_operator.precompute(x=x)

            for operator in self.residuals:
                with Builder() as m:
                    x = m.input('x')
                    wn, s, fs = replay(x)
                    r = operator(wn, s, fs)
                    m.output(y=r)

                y, [Dv1] = m.compute_with_gnDp(vout='y',
                            init=dict(x=x),
                            v=dict(x_=v))
                Dv = Dv + Dv1
            # H is 2 JtJ, see wikipedia on Gauss Newton.
            return Dv * 2

        def addmul(a, b, c, p=1):
            """ a + b * c ** p, follow the type of b """
            if p is not 1: c = c ** p
            c = b * c
            if a is not 0: c = c + a
            return c

        def dot(a, b):
            """ einsum('i,i->', a, b) """
            return self.comm.allreduce((a * b).sum())

        vs = VectorSpace(addmul=addmul, dot=dot)
        BaseProblem.__init__(self,
                        vs = vs,
                        objective=objective,
                        gradient=gradient,
                        hessian_vector_product=hessian_vector_product)

    def save(self, filename, state):
        with Builder() as m:
            x = m.input('x')
            wn, s, fs = self.forward_operator(x)
            m.output(wn=wn, s=s, fs=fs)

        wn, s, fs = m.compute(['wn', 's', 'fs'], init=dict(x=state['x']))
       
        numpy.save(filename+".wn",wn)
        numpy.save(filename+".s",s)
        numpy.save(filename+".fs",fs)

    def state(self, state):
        with Builder() as m:
            x = m.input('x')
            wn, s, fs = self.forward_operator(x)
            m.output(wn=wn, s=s, fs=fs)
        wn, s, fs = m.compute(['wn', 's', 'fs'], init=dict(x=state['x']))

        return s, fs
