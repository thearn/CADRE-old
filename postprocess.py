import pprint

import sqlite3
from cPickle import dumps, loads, HIGHEST_PROTOCOL, UnpicklingError, load, dump
import numpy as np 

from matplotlib import pylab



def get_data(db_file_name):

    conn = sqlite3.connect(db_file_name)
    conn.text_factory = sqlite3.OptimizedUnicode
    cur = conn.cursor()


    sql_count = "SELECT count(*) from cases"
    print "total cases:",  [x for x in cur.execute(sql_count)][0][0]

    sql_count = 'SELECT case_id from casevars where name=="Objective_0"'
    ids = [x[0] for x in cur.execute(sql_count)]

    id_index_map = dict([(id,i) for i,id in enumerate(ids)])
    n_cases = len(ids)
    print "top level cases:", n_cases
    data = {}

    sql = "SELECT var_id,name,case_id,sense,value FROM cases INNER JOIN casevars ON casevars.case_id=cases.id WHERE case_id in (%s)"%','.join(map(str,ids))

    cur = conn.cursor()
    cur.execute(sql)

    objective_cases = []
    for var_id, vname, case_id, sense, value in cur: 
        if not isinstance(value, (float, int, str)):
            try:
                value = loads(str(value))
            except UnpicklingError as err:
                raise UnpicklingError("can't unpickle value '%s' for"
                                      " case '%s' from database: %s"
                                      % (vname, case_id, str(err)))

        index = id_index_map[case_id]
        if vname in data: 
            data[vname][index] = value
        else: 
            data[vname]=[None,]*n_cases
            data[vname][index] = value 

    return n_cases,data


if __name__ == "__main__": 

    vars_name = {
        "Data": 1e3,
        "Gamma": np.pi/180., 
        "SOC": 1e-2, 
        "P_comm": 1/5., 
        'Comm_BitRate.gain': 1., 
        "CommLOS": 1., 
        "LOS": 1.,    
    }

    n_cases, data = get_data('CADRE.db')

    #pprint.pprint(data.keys()); exit()

    time = data['pt0.t'][0]/3600.


    iters = [0,50,150]

    for iter_index in iters: 
        for mp_index in xrange(6):
            for var, s in vars_name.iteritems(): 
                d = data['pt%d.%s'%(mp_index,var)][iter_index]
                #print var, d/s

                output = np.empty((2,time.shape[0]))
                output[0,:] = time
                output[1,:] = d/s

                np.savetxt('postproc_results/%s_%d-%d.dat'%(var,iter_index,mp_index), output.T)


#######################

    # pylab.figure()
    # pylab.plot(seconds_per_case, label="individual")
    # pylab.plot(avg_seconds_per_case, label="average")
    # pylab.legend(loc="best")
    # pylab.title('Case Computational Cost')
    # pylab.xlabel('Iteration #')
    # pylab.ylabel('Time (sec)')

    # Z = np.array(Z)
    # if not len(Z):
    #     print "no data yet..."
    #     quit()

    pylab.rcParams['text.latex.preamble']=[r"\usepackage{lmodern}"]
    params = {
        'font.size': 15, 
        'font': 'serif', 
        'font.family': "lmodern",
        'text.latex.unicode': True,
    }
    pylab.rcParams.update(params)

    sum_data = np.sum(np.array([np.array(data['pt%d.Data'%i])[:,0,-1] for i in xrange(6)]), axis=0)

    fig = pylab.figure()
    pylab.subplot(211)
    pylab.title("total data")
    pylab.plot(sum_data, 'b')
    pylab.plot([0, len(sum_data)], [3e4, 3e4], 'k--', marker="o")
    pylab.ylabel('Gigabit/sec')
    pylab.ticklabel_format(style='sci', axis='y', scilimits=(0,0))


    constraints = ["Constraint ( pt%(i)d.ConCh<=0 )", 
                   "Constraint ( pt%(i)d.ConDs<=0 )", 
                   "Constraint ( pt%(i)d.ConS0<=0 )",
                   "Constraint ( pt%(i)d.ConS1<=0 )",
                   "Constraint ( pt%(i)d.SOC[0][0]=pt%(i)d.SOC[0][-1] )"
                   ]

    c_data = []
    for con in constraints: 
        d = np.array([np.array(data[con%{'i':i}])[:,0] for i in xrange(6)])
        d = np.sum(np.ma.masked_outside(d,0,10000)**2,axis=0)**.5
        c_data.append(d)


    pylab.subplot(212)
    pylab.title(r"$\left|\left|Constraints\right|\right|_2$")
    #pylab.semilogy([0, n_cases], [0, 0], 'k--')
    #pylab.semilogy(c_data[0], marker="o", label=r"$I_{bat} - 5 \leq 0$")
    pylab.semilogy(c_data[1], label=r"$-10 - I_{bat} \leq 0$", c='c')
    pylab.semilogy(c_data[2], label=r"$0.2 - SOC \leq 0$")
    #pylab.semilogy(c_data[3], marker="o", label=r"$SOC - 1 \leq 0$")
    pylab.semilogy(c_data[4], label=r"$fSOC - iSOC = 0$")
    pylab.legend(loc="best")
    #pylab.legend(bbox_to_anchor=[1.11,1.23] ,loc="upper right")

    pylab.xlabel("Iteration #")

    fig.set_size_inches(10,7)
    fig.tight_layout()
    pylab.savefig('cadre_opt_progress.pdf', dpi=1000, bbox_inches="tight")

