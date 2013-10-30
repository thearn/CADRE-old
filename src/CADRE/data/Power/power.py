from __future__ import division
import numpy
import PowerLib
import MBI

Iph0 = 0.429
A0 = 2.66e-3
I0 = 3.86e-5
q = 1.6021e-19
n = 1.35
k = 1.3807e-23
Rs = 0.008
Rsh = 40

T1 = 50
T2 = 400
A1 = -A0*0.5 #3
A2 = A0*1.5 #10
I1 = -Iph0*1e-3
I2 = Iph0

nT = 20
nA = 20
nI = 200

T = numpy.linspace(T1, T2, nT)
A = numpy.linspace(A1, A2, nA)
I = numpy.linspace(I1, I2, nI)
Vdat = PowerLib.computevcurve(nT, nA, nI, T, A, I)

dat = numpy.zeros(3 + nT + nA + nI + nT*nA*nI)
dat[:3] = [nT, nA, nI]
dat[3:3+nT] = T[:]
dat[3+nT:3+nT+nA] = A[:]
dat[3+nT+nA:3+nT+nA+nI] = I[:]
dat[3+nT+nA+nI:] = Vdat.reshape((nT*nA*nI),order='F')

numpy.savetxt('curve.dat', dat)
