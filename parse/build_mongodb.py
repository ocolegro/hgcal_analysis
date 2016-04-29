import sys
import numpy as np
from pymongo import MongoClient
sys.path.insert(0, '../analyze')
import utils as u
from multiprocessing import Pool
import time





reset_db,f_events,n_cores,n_runs = True,2500,6,10

pool     = Pool(processes=n_cores)              # start 4 worker processes


p_merged,merged_energy,p_unmerged,unmerged_energy  \
                ='gamma_','30','electron_','4000'


runs            = np.array(range(1,n_runs)).astype(str);
print 'Loading Merged events now'

merged          = u.load_events(merged_energy,p_merged,runs,f_events,'merged_txts/',n_cores)
print 'Inserting Merged events'
client   = MongoClient()
db       = client['hgcal']
if reset_db:
    db.merged.remove()

for key in merged.keys():
    db.merged.insert(merged[key])


print 'Loading Unmerged events now'
runs            = np.array(range(1,n_runs)).astype(str)
unmerged        = u.load_events(unmerged_energy,p_unmerged,runs,f_events,'unmerged_txts/',n_cores)

print 'Inserting Unmerged events'
if reset_db:
    db.unmerged.remove()

for key in unmerged.keys():
    db.unmerged.insert(unmerged[key])
