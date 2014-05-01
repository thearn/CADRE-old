import time

from openmdao.lib.casehandlers.api import CSVCaseRecorder, DBCaseRecorder
from CADRE import CADRE_Optimization

print "setting up"
top = CADRE_Optimization(n=1500, m=300)
#top.recorders = [CSVCaseRecorder(filename='CADRE.csv')]

top.recorders = [DBCaseRecorder(dbfile='CADRE.db', append=False)]
printvars = []
for var in ['Data', 'ConCh', 'ConDs', 'ConS0', 'ConS1', 'SOC']:
   printvars += ["pt" + str(i) + ".Data" for i in xrange(6)]

top.printvars = printvars
print "running"

start_time = time.time()

top.run()

print "Finished"
print "Run Time: ", time.time() - start_time