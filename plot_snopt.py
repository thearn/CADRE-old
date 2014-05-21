import numpy as np
import pylab as p

f = open('SNOPT_print.out')

data =[]

next_line = False
for line in f: 
    if "Itns Major Minors    Step   nCon Feasible  Optimal  MeritFunction     L+U BSwap     nS  condHz Penalty" in line: 
        next_line = True
    elif next_line: 
        row = line.strip().split()
        try: 
            int(row[4])
            data.append(map(float,(row[5], row[6], row[7])))
        except: 
            data.append(map(float,(row[4], row[5], row[6])))
        next_line=False



data = np.array(data)

p.plot(np.log10(data[:,0]), label="feasibility")
p.plot(np.log10(data[:,1]), label="optimality")
p.legend(loc="best")
p.title('Optimization Progress')

p.figure()
p.plot(data[:,2], label="merit function")
p.title('Objective Function')

p.show()

