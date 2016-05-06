import pickle as pkl
import numpy  as np
import matplotlib.pyplot as plt
import utils as u
from pymongo import MongoClient
import random

def merged_eng(eng_1,eng_2,dr):
    return eng_1


def mask_eng(eng_1,eng_2,dr):
    return eng_1


def smear_eng(eng,pdgid):
    return eng


def propogate(x,y,phi,theta,distance):
    r = distance* np.tan(theta)
    x_ret = x + r*np.cos(phi)
    y_ret = y + r*np.sin(phi)

    return x_ret,y_ret
'''
def propogate(x,y,phi,theta,distance):
    return x,y
'''
def sep_matrix(tmp_):
    ret_mat = []
    for particle in tmp_:
        x = particle[0]
        y = particle[1]
        sep = np.sqrt(np.power((tmp_[:,0] - x),2) + np.power((tmp_[:,1] - y),2))
        sep[sep == 0] =  1e6
        ret_mat.append(sep)
    return np.array(ret_mat)


def build_events(beam_std,bunch_size,n_events,detector_distance,cut,dead_time=5):
    client      = MongoClient()
    db          = client['hgcal']
    events      = []
    unmerg,merg,merg_eng,phot_merg_eng = [],[],[],[]
    for counter,item in enumerate(db.target.find().limit( n_events )):
        tmp_ = []
        for id,particle in enumerate(item['particles']):
            if particle['eng'] < cut : continue

            #x_vtx,y_vtx   = particle['x_pos']+np.random.normal(0,beam_std),particle['y_pos']+np.random.normal(0,beam_std)
            x_vtx,y_vtx   = np.random.normal(0,beam_std),np.random.normal(0,beam_std)
            x_cord,y_cord = propogate(x_vtx,y_vtx,particle['phi'],particle['theta'],detector_distance)
            tmp_.append([x_cord,y_cord,particle['pdg'],particle['eng'],dead_time])
        if len(tmp_) == 0:  continue

        if len(events) > 0: events = np.vstack((tmp_,events[events[:,-1] > 0]))
        else: events = np.array(tmp_)

        if counter < 10: continue
        if counter%bunch_size == 0:
            if len(events) > 0: events[:,-1] = events[:,-1]  - 1

            #sep_mat = sep_matrix(events)
            evt_sub = events[events[:,-1]==dead_time-1.0]
            sep_sub = sep_matrix(evt_sub)
            for id,ele in enumerate(evt_sub):

                if evt_sub[id][2] != 22.: continue
                argmin = np.argmin(sep_sub[id])
                if sep_sub[id][argmin] > 20.:
                    unmerg.append(1)
                else:
                    merg.append(1)
                    merg_eng.append(evt_sub[argmin][3] )
                    phot_merg_eng.append(evt_sub[id][3])
                #dr.append(sep_sub[id][argmin])

    return np.array(merg),np.array(unmerg)
detector_distance,beam_std,bunch_size,measurement_cut,n_events\
               = 250,21,2,2000,2000

sim            = pkl.load( open( "../output/pkl/eff.pkl", "rb" ) )
in_eff         = 1-sim['eff_tuple'][2]
eff = []
for detector_distance in [2.5,5.,7.5,10.,12.5,15.,17.5,20.,22.5,25.,27.5,30.,32.5,35.]:
    merg,unmerg       = build_events(beam_std,bunch_size,n_events,detector_distance,measurement_cut)
    eff.append(float(len(merg))/(len(merg)+len(unmerg)))

plt.plot([2.5,5.,7.5,10.,12.5,15.,17.5,20.,22.5,25.,27.5,30.,32.5,35.],eff)
plt.xlabel('Detector Distance',fontsize=20)
plt.ylabel('Fraction of Merged Photons',fontsize=20)
plt.savefig('../output/plots/merge_vs_std.png')
plt.show()