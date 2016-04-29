import pickle as pkl
import numpy  as np
import matplotlib.pyplot as plt


def identity(x):
    return x


def binned_vals(bins_new,binned_dist,bin_multiplier,trans_func = identity,spacing = 0):
    ret_dist,binned_dist = [],np.array(binned_dist)

    for id,bin in enumerate(bins_new):
        if id           == len(bins_new) -1 : break
        bin_low,bin_high = bin,bins_new[id+1]

        temp_bin        = binned_dist[[id for id,x in enumerate(abs(spacing - trans_func(binned_dist[:,0]) * bin_multiplier)) if bin_low < x  < bin_high ]]
        ret_dist.append(np.sum(temp_bin[:,1]))
    return ret_dist


gev_to_mev,detector_distance,electron_spacing,neighbor_prob\
               = 1000,350,25,.01
'''
The goal is to simulate the probability that a bremm'ed photon overlaps with a nearby electron and is undetected.


(1) gev_to_mev is a simple conversion to help relate binning
(2) detector_distance is the proposed distance from the detector to the target, in mm.
(3) electron_spacing is the spacing between the incoming electron and it's nearby neighbors
(4) neighbor_prob  is the probability for a random point picked at a radius of electron_spacing to contain a nearby beam-line electron/

Radial prob is something I cooked up because even if a hard photon radiates at the distance of electron spacing from the inc. electron it will still need to overlap with a nearby neighbor.

'''
sim            = pkl.load( open( "sim.pkl", "rb" ) )
eng_dist       = binned_vals(sim['eff_tuple'][1],sim['eng_dist'],gev_to_mev)
print np.array(sim['ang_dist'])
sep_dist       = binned_vals(sim['eff_tuple'][0],sim['ang_dist'],detector_distance,np.sin,electron_spacing)
missed_frac    = neighbor_prob*np.dot(sep_dist,np.dot(eng_dist,1-sim['eff_tuple'][2]))
print  'The fraction of photons which merge undetected is %s ' %(missed_frac)

bins = np.linspace(1,30,90)

eff_matrix = np.dot(eng_dist,1-sim['eff_tuple'][2])

plt.semilogy(bins,[np.dot(binned_vals(sim['eff_tuple'][0],sim['ang_dist'],detector_distance,np.sin,distance),eff_matrix) for distance in bins],'ro')
plt.xlabel('Electron Separation, r',size=20)
plt.ylabel('Merging Probabiilty/Electron Density',size=20)
plt.title('Geant4 Sim.',size=16,loc='left')
plt.savefig('ele_rejection.png')