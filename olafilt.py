#
# Overlap-add FIR filter, (c) Joachim Thiemann 2016
# LICENCE: MIT
#
import numpy as np

def olafilt(b, x, zi=None):
    """
    Filter a one-dimensional array with an FIR filter

    Filter a data sequence, `x`, using a FIR filter given in `b`.
    Filtering uses the overlap-add method converting both `x` and `b`
    into frequency domain first.  The FFT size is determined as the
    next higher power of 2 of twice the length of `b`.

    Parameters
    ----------
    b : one-dimensional numpy array
        The impulse response of the filter
    x : one-dimensional numpy array
        Signal to be filtered
    zi : one-dimensional numpy array, optional
        Initial condition of the filter, but in reality just the
        runout of the previous computation.  If `zi` is None or not
        given, then zero initial state is assumed.

    Returns
    -------
    y : array
        The output of the digital filter.
    zf : array, optional
        If `zi` is None, this is not returned, otherwise, `zf` holds the
        final filter delay values.
    """

    L_I = b.shape[0]
    # Find power of 2 larger that 2*L_I (from abarnert on Stackoverflow)
    L_F = 2<<(L_I-1).bit_length()  
    L_S = L_F - L_I + 1
    L_sig = x.shape[0]
    offsets = range(0, L_sig, L_S)

    # blockwise frequency domain multiplication
    if np.iscomplexobj(b) or np.iscomplexobj(x):
        FDir = np.fft.fft(b, n=L_F)
        tempresult = [np.fft.ifft(np.fft.fft(x[n:n+L_S], n=L_F)*FDir)
                      for n in offsets]
        res = np.zeros(L_sig+L_F, dtype=np.complex128)
    else:
        FDir = np.fft.rfft(b, n=L_F)
        tempresult = [np.fft.irfft(np.fft.rfft(x[n:n+L_S], n=L_F)*FDir)
                      for n in offsets]
        res = np.zeros(L_sig+L_F)

    # overlap and add
    for i, n in enumerate(offsets):
        res[n:n+L_F] += tempresult[i]

    if zi is not None:
        res[:zi.shape[0]] = res[:zi.shape[0]] + zi
        return res[:L_sig], res[L_sig:]
    else:
        return res[:L_sig]
