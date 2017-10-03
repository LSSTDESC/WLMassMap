# This module contains the code for a simple flat Kaiser-Squires inversion
import numpy as np


def flat_KS_map(gmap):
    """Compute kappa maps from binned shear maps.
    returns kappa_e and kappa_b
    """
    nx = gmap.shape[1]
    ny = gmap.shape[2]

    # get frequencies
    k1, k2 = np.meshgrid(np.fft.fftfreq(nx),
                         np.fft.fftfreq(ny))

    g1 = np.fft.fft2(gmap[0])
    g2 = np.fft.fft2(gmap[1])
    denom = k1*k1 + k2*k2
    denom[0, 0] = 1  # avoid division by 0
    kap = ((k1*k1 - k2*k2) - 2j*(k1*k2)) * (g1 + 1j*g2) / denom
    kap = np.fft.ifft2(kap)

    return np.real(kap), np.imag(kap)
