from __future__ import division
import numpy, pylab
import PowerLib
import MBI

Iph0 = 0.429
A0 = 1.834e-3
I0 = 3.86e-5
q = 1.6021e-19
n = 1.35
k = 1.3807e-23
Rs = 0.008
Rsh = 40

T1 = 50
T2 = 400
A1 = 0
A2 = A0

V1 = -25
V2 = 25
I1 = 0
I2 = Iph0

nT = 7
nA = 2
n = 5000
T = numpy.linspace(T1, T2, nT)
A = numpy.linspace(A1, A2, nA)
V = numpy.linspace(V1, V2, n)
I = numpy.linspace(I1, I2, n)

Idat = PowerLib.computeivcurve(nT, nA, n, Iph0, A0, I0, q, n, k, Rs, Rsh, T, A, V)
Vdat = numpy.zeros((nT, nA, n), order='F')
Ieval = numpy.zeros((n,1),order='F')

for iA in range(nA):
    for iT in range(nT):
        m = MBI.MBI(V[:], [Idat[iT,iA,:]], [int(n/3)], [4])
        maxI = max(Idat[iT,iA,:])
        minI = min(Idat[iT,iA,:])
        Ieval[:,0] = I[:]
        for i in range(n):
            if Ieval[i,0] > maxI:
                Ieval[i,0] = maxI
            elif Ieval[i,0] < minI:
                Ieval[i,0] = minI
        Vdat[iT,iA,:] = m.evaluate(Ieval)[:,0]
        print 'T:', iT, 'A:', iA, 'I range:', min(Idat[iT,iA,:]), max(Idat[iT,iA,:])

Pdat = numpy.zeros((nA, nA, nA, nA, nA, nA, nA, nT), order='F')
for iT in range(nT):
    Pdat[:,:,:,:,:,:,:,iT] = PowerLib.computemaxpower(nA, n, I, Vdat[iT,:,:])

numpy.savetxt('power.dat',Pdat.flatten(order='F'))

info = numpy.zeros(2 + nT + nA)
info[0] = nT
info[1] = nA
info[2:2+nT] = T[:]
info[2+nT:2+nT+nA] = A[:]

numpy.savetxt('powerInfo.dat', info)
