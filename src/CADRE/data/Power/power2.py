from __future__ import division
import numpy
import PowerLib
import MBI
import scipy.optimize

V0 = -0.6
AT = 2.66e-3
Isc0 = 0.453
q = 1.6021e-19
n = 1.35
k = 1.3807e-23
Rs = 0.008
Rsh = 40
Isat = 2.809e-12

T1 = 50
T2 = 400
A1 = -AT*0.5 #3
A2 = AT*1.5 #10
I1 = -Isc0*1e-3
I2 = Isc0

nT = 20
nA = 20
nI = 200

T = numpy.linspace(T1, T2, nT)
A = numpy.linspace(A1, A2, nA)
I = numpy.linspace(I1, I2, nI)

def M(Isc,VT,I):
    f = lambda V: Isc + Isat - I - Isat*numpy.exp(V/VT) - V/Rsh
    if I < Isc:
        V = scipy.optimize.brentq(f, -5, 5)
    else:
        V = V0 * numpy.tanh(-VT*Rsh/V0/(Isat*Rsh+VT)*(I-Isc))
    return V
#Vdat = PowerLib.computevcurve(nT, nA, nI, T, A, I)

Vdat = numpy.zeros((nT,nA,nI),order='F')
for iT in range(nT):
    for iA in range(nA):
        for iI in range(nI):
            Isc = A[iA]/AT*Isc0
            VT = n*k*T[iT]/q
            Vdat[iT,iA,iI] = M(Isc,VT,I[iI])

dat = numpy.zeros(3 + nT + nA + nI + nT*nA*nI)
dat[:3] = [nT, nA, nI]
dat[3:3+nT] = T[:]
dat[3+nT:3+nT+nA] = A[:]
dat[3+nT+nA:3+nT+nA+nI] = I[:]
dat[3+nT+nA+nI:] = Vdat.reshape((nT*nA*nI),order='F')

numpy.savetxt('curve.dat', dat)
