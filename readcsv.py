import csv
import numpy
 as np
import pylab

f = open("CADRE.csv", "rb")
reader = csv.DictReader(f, skipinitialspace=True)

X, Y, Z = [], [], []

for row in reader:
    data = [row["pt" + str(i) + ".Data[0][1499]"] for i in xrange(6)]
    sumdata = sum([float(i) for i in data if i])
    c1 = [row["Constraint ( pt" + str(i) + ".ConCh<=0 )"] for i in xrange(6)]
    c2 = [row["Constraint ( pt" + str(i) + ".ConDs<=0 )"] for i in xrange(6)]
    c3 = [row["Constraint ( pt" + str(i) + ".ConS0<=0 )"] for i in xrange(6)]
    c4 = [row["Constraint ( pt" + str(i) + ".ConS1<=0 )"] for i in xrange(6)]
    c5 = [row["Constraint ( pt" + str(i) + ".SOC[0][0]=pt" + str(i) + ".SOC[0][-1] )"]
          for i in xrange(6)]
    # c1_f = np.all([float(i) < 0 for i in c1 if i])
    # c2_f = np.all([float(i) < 0 for i in c2 if i])
    # c3_f = np.all([float(i) < 0 for i in c3 if i])
    # c4_f = np.all([float(i) < 0 for i in c4 if i])
    # c5_f = np.all([float(i) < 0 for i in c4 if i])

    c1_f = sum([float(i) for i in c1 if i])
    c2_f = sum([float(i) for i in c2 if i])
    c3_f = sum([float(i) for i in c3 if i])
    c4_f = sum([float(i) for i in c4 if i])
    c5_f = sum([float(i) for i in c5 if i])

    feasible = [c1_f, c2_f,  c3_f, c4_f, c5_f]

    X.append(sumdata), Y.append(sum(feasible)), Z.append(feasible)

    # print sumdata, sum(feasible), max(feasible) #,[ '%.1f' % i for i in
    # feasible]
    print sumdata


Z = np.array(Z)


pylab.subplot(311)
pylab.title("total data")
pylab.plot(X, 'b')
pylab.plot([0, len(X)], [3e4, 3e4], 'k--')
pylab.subplot(312)
pylab.title("Sum of Constraints")
pylab.plot([0, len(Y)], [0, 0], 'k--')
pylab.plot(Y, 'k')
pylab.subplot(313)
pylab.title("Max of Constraints")
pylab.plot([0, len(Z)], [0, 0], 'k--')
pylab.plot(Z[:, 0], label="c1")
pylab.plot(Z[:, 1], label="c2")
pylab.plot(Z[:, 2], label="c3")
pylab.plot(Z[:, 3], label="c4")
pylab.plot(Z[:, 4], label="c5")
pylab.legend(loc="best")
pylab.show()
