import sys
import numpy as np
from pymongo import MongoClient
sys.path.insert(0, '../analyze')
import utils as u
from multiprocessing import Pool
import time





reset_db,f_events,n_cores,n_runs,step = False,2500,6,1000,24
single_gamma,single_ele,ele_gamma,ele_ele                = False,True,False,False
client   = MongoClient()
db       = client['hgcal']
pool     = Pool(processes=n_cores)              # start 4 worker processes

if single_gamma:
    if reset_db:
        db.photons.remove()
    for prefix in ['gamma_single_30_']:
        for i in range(int(n_runs/step)):

            runs            = np.array(range(1 + step*i,step*(i+1))).astype(str);
            print 'Loading Merged events now'

            merged          = u.load_events('unmerged_files/',prefix,runs,f_events,n_cores)

            print 'Inserting Merged events'



            for key in merged.keys():
                try:
                    if set(merged[key]['gen_pdgid']) != set([22.]):
                        continue
                    db.photons.insert(merged[key])
                except: pass

if single_ele:
    if reset_db:
        db.electrons_finer.remove()
    for prefix in ['electron_single_30_']:
        for i in range(int(n_runs/step)):

            runs            = np.array(range(1 + step*i,step*(i+1))).astype(str);
            print 'Loading Merged events now'

            merged          = u.load_events('unmerged_files/',prefix,runs,f_events,n_cores)

            print 'Inserting Merged events'



            for key in merged.keys():
                try:
                    if set(merged[key]['gen_pdgid']) != set([11.]):
                        continue
                    db.electrons_finer.insert(merged[key])
                except: pass

if ele_gamma:

    if reset_db:
        db.ele_gamma.remove()
    for prefix in ['electron_photon_30_30_','electron_single_30_','photon_single_30_','electron_electron_30_']:
        for i in range(int(n_runs/step)):
            print prefix
            runs            = np.array(range(1 + step*i,step*(i+1))).astype(str);
            print 'Loading Merged events now'

            merged          = u.load_events('merged_files/',prefix,runs,f_events,n_cores)

            print 'Inserting Merged events'


            for key in merged.keys():
                try:
                    if set(merged[key]['gen_pdgid']) != set([11.,22.]):
                        continue
                    db.ele_gamma_finer.insert(merged[key])
                except: pass

if ele_ele:

    if reset_db:
        db.ele_ele.remove()
    for prefix in ['electron_electron_30_','electron_single_30_']:
        for i in range(int(n_runs/step)):
            print prefix
            runs            = np.array(range(1 + step*i,step*(i+1))).astype(str);
            print 'Loading Merged events now'

            merged          = u.load_events('merged_files/',prefix,runs,f_events,n_cores)

            print 'Inserting Merged events'


            for key in merged.keys():
                try:
                    if set(merged[key]['gen_pdgid']) != set([11.,11.]):
                        continue
                    db.ele_ele.insert(merged[key])
                except: pass


'''
    print 'Loading Unmerged events now'
    runs            = np.array(range(1 + step*i,step*(i+1))).astype(str)
    unmerged        = u.load_events(unmerged_energy,p_unmerged,runs,f_events,'unmerged_txts_2/',n_cores)

    print 'Inserting Unmerged events'
    if reset_db:
        db.unmerged.remove()

    for key in unmerged.keys():
        db.unmerged.insert(unmerged[key])
'''
