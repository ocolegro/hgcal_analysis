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
    gamma_x,gamma_y,gamma_eng,ele_x,ele_y,ele_eng,merg_x,merg_y,merg_eng,merg_peng,merg_eeng,merg_dr,unmerg_counter = [],[],[],[],[],[],[],[],[],[],[],[],[]
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

            evt_sub = events[events[:,-1]==dead_time-1.0]
            sep_sub = sep_matrix(evt_sub)
            for id,ele in enumerate(evt_sub):
                if evt_sub[id][2] == 22.:
                    gamma_x.append(evt_sub[id][0])
                    gamma_y.append(evt_sub[id][1])
                    gamma_eng.append(evt_sub[id][3])
                    ##Check for merging

                    argmin = np.argmin(sep_sub[id])
                    #    unmerg.append(1)
                    #else:
                    #    merg.append(1)
                    dr = sep_sub[id][argmin]
                    if  dr < 20.:
                        #For now, we ignore photon-photon mergers
                        if evt_sub[argmin][2] == 22.: continue
                        merg_x.append(evt_sub[id][0] )
                        merg_y.append(evt_sub[id][1])
                        merg_eng.append(evt_sub[argmin][3] + evt_sub[id][3])
                        merg_peng.append(evt_sub[id][3])
                        merg_eeng.append(evt_sub[argmin][3])

                        merg_dr.append(dr)
                    else:
                        unmerg_counter.append(1)

                if evt_sub[id][2] == 11.:
                    ele_x.append(evt_sub[id][0])
                    ele_y.append(evt_sub[id][1])
                    ele_eng.append(evt_sub[id][3])


    return np.array(gamma_x),np.array(gamma_y),np.array(gamma_eng),np.array(ele_x),np.array(ele_y),np.array(ele_eng),np.array(merg_x),np.array(merg_y),np.array(merg_eng),\
           np.array(merg_peng),np.array(merg_eeng),np.array(merg_dr),np.array(unmerg_counter)
detector_distance,beam_std,bunch_size,measurement_cut,n_events\
               = 250,21,2,2000,2000000

sim            = pkl.load( open( "../output/pkl/eff.pkl", "rb" ) )
in_eff         = 1-sim['eff_tuple'][2]
gamma_x,gamma_y,gamma_eng,ele_x,ele_y,ele_eng,merg_x,merg_y,merg_eng,\
merg_peng,merg_eeng,merg_dr,unmerg_counter       \
               = build_events(beam_std,bunch_size,n_events,detector_distance,measurement_cut)


plt.hist(ele_eng)
plt.yscale('log', nonposy='clip',bins=30)
plt.xlabel('Transmitted Electron Energy',fontsize=20)
plt.ylabel('Nevents',fontsize=20)
plt.title('Geant4 Sim. (%s entries)' %(len(ele_eng)),fontsize=12,loc='left')
plt.savefig('../output/plots/trans_ele_eng.png')
plt.clf()


plt.hist(gamma_eng,color='red')
plt.xlabel('Transmitted Photon Energy',fontsize=20)
plt.ylabel('Nevents',fontsize=20)
plt.title('Geant4 Sim. (%s entries)' %(len(gamma_eng)),fontsize=12,loc='left')
plt.savefig('../output/plots/trans_gamma_eng.png')
plt.clf()

plt.hist(merg_eng,color='red')
plt.xlabel('Total Merged Energy',fontsize=20)
plt.ylabel('Nevents',fontsize=20)
plt.title('Geant4 Sim. (%s entries)' %(len(merg_eng)),fontsize=12,loc='left')
plt.savefig('../output/plots/merg_eng.png')
plt.clf()

plt.hist(merg_eng,color='maroon')
plt.xlabel('Total Merged Energy',fontsize=20)
plt.ylabel('Nevents',fontsize=20)
plt.title('Geant4 Sim. (%s entries)' %(len(merg_eng)),fontsize=12,loc='left')
plt.savefig('../output/plots/merg_eng.png')
plt.clf()

plt.hist(merg_eeng,color='blue')
plt.yscale('log', nonposy='clip',bins=30)
plt.xlabel('Merged Electron Energy',fontsize=20)
plt.ylabel('Nevents',fontsize=20)
plt.title('Geant4 Sim. (%s entries)' %(len(merg_eeng)),fontsize=12,loc='left')
plt.savefig('../output/plots/merg_eeng.png')
plt.clf()

plt.hist(merg_peng,color='red')
plt.xlabel('Merged Photon Energy',fontsize=20)
plt.ylabel('Nevents',fontsize=20)
plt.title('Geant4 Sim. (%s entries)' %(len(merg_peng)),fontsize=12,loc='left')
plt.savefig('../output/plots/merg_peng.png')
plt.clf()



plt.hist(merg_dr,color='maroon')
plt.xlabel('Merged Photon $\\Delta R$',fontsize=20)
plt.ylabel('Nevents',fontsize=20)
plt.title('Geant4 Sim. (%s entries)' %(len(merg_dr)),fontsize=12,loc='left')
plt.savefig('../output/plots/merg_dr.png')
plt.clf()


plt.scatter(ele_x[0:len(gamma_x)],ele_y[0:len(gamma_y)],marker ='.')
plt.xlabel('X-Position',fontsize=20)
plt.ylabel('Y-Position',fontsize=20)
plt.title('Geant 4 Sim: Electron POI (%s entries)' %(len(ele_x)),fontsize=12,loc='left')
plt.xlim(-80,80)
plt.ylim(-80,80)
plt.savefig('../output/plots/ele_scatter.png')
plt.clf()


plt.scatter(gamma_x,gamma_y,marker ='.',color='red')
plt.xlabel('X-Position',fontsize=20)
plt.ylabel('Y-Position',fontsize=20)
plt.title('Geant 4 Sim: Photon POI (%s entries)' %(len(gamma_x)),fontsize=12,loc='left')
plt.xlim(-80,80)
plt.ylim(-80,80)
plt.savefig('../output/plots/gamma_scatter.png')
plt.clf()

plt.scatter(ele_x,ele_y,marker ='.')
plt.scatter(gamma_x,gamma_y,marker ='.',color='red')
plt.xlabel('X-Position',fontsize=20)
plt.ylabel('Y-Position',fontsize=20)
plt.title('Geant 4 Sim: Superposed Electron-Photon POI (%s entries)' %(len(ele_x)+len(gamma_x)),fontsize=12,loc='left')
plt.xlim(-80,80)
plt.ylim(-80,80)
plt.savefig('../output/plots/super_scatter.png')
plt.clf()


plt.scatter(merg_x,merg_y,marker ='.',color='maroon')
plt.xlabel('X-Position',fontsize=20)
plt.ylabel('Y-Position',fontsize=20)
plt.title('Geant 4 Sim: Merged POI (%s entries)' %(len(merg_x)),fontsize=12,loc='left')
plt.xlim(-80,80)
plt.ylim(-80,80)
plt.savefig('../output/plots/merg_scatter.png')
plt.clf()



'''

eff = []
for detector_distance in np.arange(100,500,50):
    gamma_x,gamma_y,gamma_eng,ele_x,ele_y,ele_eng,merg_x,merg_y,merg_eng,merg_dr,unmerg_counter       = build_events(beam_std,bunch_size,n_events,detector_distance,measurement_cut)
    eff.append(float(len(merg_dr))/(len(merg_dr)+len(unmerg_counter)))

plt.plot([2.5,5.,7.5,10.,12.5,15.,17.5,20.,22.5,25.,27.5,30.,32.5,35.],eff)
plt.xlabel('Detector Distance',fontsize=20)
plt.ylabel('Fraction of Merged Photons',fontsize=20)
plt.savefig('../output/plots/merge_vs_dist.png')
plt.show()


eff = []
for beam_std in np.arange(5,50,5):
    gamma_x,gamma_y,gamma_eng,ele_x,ele_y,ele_eng,merg_x,merg_y,merg_eng,merg_dr,unmerg_counter       = build_events(beam_std,bunch_size,n_events,detector_distance,measurement_cut)
    eff.append(float(len(merg_dr))/(len(merg_dr)+len(unmerg_counter)))

plt.plot([2.5,5.,7.5,10.,12.5,15.,17.5,20.,22.5,25.,27.5,30.,32.5,35.],eff)
plt.xlabel('Detector Distance',fontsize=20)
plt.ylabel('Fraction of Merged Photons',fontsize=20)
plt.savefig('../output/plots/merge_vs_std.png')
plt.show()
'''