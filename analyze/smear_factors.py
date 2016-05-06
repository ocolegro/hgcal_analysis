import sys
sys.path.insert(0,"../")
import utils as u
import smear_functions as sf
from pymongo import MongoClient
import matplotlib.pyplot as plt

n_events  = 200000
##Load up specified number of events from MongoDB
client      = MongoClient()
db          = client['hgcal']




arr       = u.prep_smear(u.load_event_array(db.electrons_finer.find(no_cursor_timeout=True).limit(n_events)),False,[11.])
sf.single_smear(arr[:,2],arr[:,3],'single_electron')

#arr       = u.prep_smear(u.load_event_array(db.photons.find(no_cursor_timeout=True).limit(n_events)),False,[22.])
#sf.single_smear(arr[:,1],arr[:,3],'single_photon')

arr = u.prep_smear(u.load_event_array(db.ele_gamma_finer.find(no_cursor_timeout=True).limit(n_events)),True,[11.,22.])
sf.merged_smear(arr[:,1],arr[:,2],arr[:,3],'electron_photon')



'''
arr = u.prep_mask(u.load_event_array(db.photons.find(no_cursor_timeout=True).limit(100000)),
                   u.load_event_array(db.electrons.find(no_cursor_timeout=True).limit(100000)),n_events)
print 'Done'
plt.scatter(arr[:,0],arr[:,2])
plt.show()
'''