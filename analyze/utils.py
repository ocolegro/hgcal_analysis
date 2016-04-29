import numpy as np
import matplotlib.pyplot as plt
import utils as u
import pickle as pkl
from sklearn.metrics import roc_curve, auc
from decimal import Decimal


def weights(x):
        if x < 8:
            y = 1.0
        elif x < 19:
            y = 19.913/15.7789
        else:
            y = 26.9505/15.7789
        return y


def load_file(arr_):
    event = arr_[2] * arr_[3]
    try:
        lines = [line.rstrip('\n') for line in open(arr_[0])]
        for line in lines:
            try:
                if line.find('>') > 0:
                    continue
                var = str(line[1:line.find("=")-1]);arr = line[line.find("=")+2:]
                if 'Event' in var:
                    continue

                if 'event_' in line:
                    event = event + 1
                    arr_[1][event] = {}
                    continue

                if '' in arr.split(","):break
                low_  = var.find(".");high_ = var.find("_")

                if var[0:low_] == 'HGCSSGenParticleVec':
                    key = 'gen_'  + var[low_+1:high_]
                else:
                    key = 'reco_' + var[low_+1:high_]
                arr_[1][event][key] = [float(i) for i in arr.split(",")]
            except:
                pass
    except:
        print 'The file %s does not exist' % (arr_[0])
    return arr_[1]


def load_events(energy,process,runs,chunk_size,suffix = '',n_cores = 6):
        from multiprocessing import Pool
        pool  = Pool(processes=n_cores)              # start 4 worker processes
        print 'Pooling inputs'
        inputs = [ ["../data/" + suffix +  process + energy + "_" + run + "_MEV.txt",{},chunk_size,int(run)] for run in runs]
        print 'Inputs are loaded, creating super_dict'
        dicts = pool.map(load_file, inputs)
        super_dict = {}
        for k in set(k for d in dicts for k in d):
            super_dict[str(k)] = [d[k] for d in dicts if k in d][0]
        return super_dict


def prep(e_,bkg = False):
    keys    = e_.keys()
    wgt     = np.vectorize(u.weights)
    ret     = []
    for key in keys:
        if len(e_[key].keys()) != 23: continue
        for key_2 in e_[key].keys(): e_[key][key_2] = np.array(e_[key][key_2])

        if len(e_[key]['gen_pz']) == 0: continue
        if bkg:
            if set(np.array(e_[key]['gen_pdgid'])) != set([11.,22.]): continue
            if e_[key]['gen_pz'][np.array(e_[key]['gen_pdgid'])==22.][0] < 1900: continue
            if e_[key]['gen_pz'][np.array(e_[key]['gen_pdgid'])==22.][0] > 4100: continue
            if np.sqrt(  np.power((e_[key]['gen_ypos'][0] -e_[key]['gen_ypos'][1]),2)) > 10: continue
        else:
            if len(e_[key]['gen_pdgid']) != 1:
                continue
            if e_[key]['gen_pdgid'] != [11.] :
                continue
            if e_[key]['gen_pz'][0] < 3900:
                continue
            if e_[key]['gen_pz'][0] > 4100: continue

        if np.max(np.abs(e_[key]['reco_xpos'])) > 100: continue
        if np.max(np.abs(e_[key]['reco_ypos'])) > 100: continue

        lay_s,eng_s = e_[key]['reco_layer'],e_[key]['reco_energy']
        len_s       = len(lay_s)

        eng_dot       = np.dot(wgt(lay_s),eng_s)
        if eng_dot != eng_dot:continue
        x_s,y_s     = e_[key]['reco_xpos'],e_[key]['reco_ypos']

        if len(x_s) != len(y_s): continue
        min_s,max_s = np.min(lay_s),np.max(lay_s)
        if min_s == 0: min_s = 1
        range_s = min_s-max_s

        frac_hit_s,frac_eng_s  = [],[]

        for ind in range(1,4):
            frac_hit_s.append(float(len(lay_s[lay_s == ind]))/len_s)
            frac_eng_s.append(np.sum(eng_s[lay_s == ind])/eng_dot)

        mol_s1,mol_s2       = [],[]
        for line in range(6,28):
            try:
                mol1,mol2 = u.calc_cont_r(x_s,y_s,eng_s,lay_s,1,line)
            except:
                mol1,mol2 = 0,0
            mol_s1.append(mol1);mol_s2.append(mol2);

        if bkg == False:
            vars_s        = [1,-10000,len_s,eng_dot,mol_s1[1],mol_s1[-1],mol_s2[1],mol_s2[-1],np.append(frac_hit_s,frac_eng_s),range_s,max_s,min_s]

        else:
            separation    = np.sqrt( np.power(e_[key]['gen_xpos'][0]-e_[key]['gen_xpos'][1],2) + np.power(e_[key]['gen_ypos'][0]-e_[key]['gen_ypos'][1],2))
            eng           = e_[key]['gen_pz'][np.array(e_[key]['gen_pdgid'])==22.][0]
            vars_s        = [separation,eng,len_s,eng_dot,mol_s1[1],mol_s1[-1],mol_s2[1],mol_s2[-1],np.append(frac_hit_s,frac_eng_s),range_s,max_s,min_s]

        vars_s = np.hstack(vars_s)
        ret.append(vars_s);

    return np.array(ret)


def load_event_array(db_merged,db_unmerged):
    merged_evts,unmerged_evts = {},{}
    for id,event in enumerate(db_merged):
        merged_evts[id] = (event)
    merged_evts = u.prep(merged_evts,True)

    for id,event in enumerate(db_unmerged):
        unmerged_evts[id] = event
    unmerged_evts = u.prep(unmerged_evts)
    return merged_evts,unmerged_evts


def calc_cont_r(x,y,e_layer,layers,cut_l=0,cut_h=100,confinement = [.68,.90]):
    cut_l = layers >= cut_l

    x = x[ cut_l];y = y[cut_l];e_layer = e_layer[cut_l];layers_l = layers[cut_l]

    cut_h = layers_l  <= cut_h

    if len(layers_l[cut_h]) == 0: return 0,0

    x = x[cut_h];y = y[cut_h];e_layer = e_layer[cut_h]
    x_mu  = np.sum(x*e_layer/np.sum(e_layer));y_mu = np.sum(y*e_layer/np.sum(e_layer))
    radii = np.sqrt((x-x_mu)*(x-x_mu) + (y - y_mu)*(y - y_mu))
    inds  = np.argsort(radii)
    e_    = e_layer[inds];radii = radii  [inds]
    e_    = np.cumsum(e_);e_ = e_/e_[-1]
    ret,perc = [],[]

    if len(x) >= 3:
        for line in confinement:
            percentile = np.searchsorted(e_,line)
            mu = radii[percentile]
            ret.append(mu)
    else:
        ret.append(np.sum(radii*e_layer/np.sum(e_layer)));ret.append(np.mean(radii*e_layer/np.sum(e_layer)))

    if ret[0] != ret[0]:
        return 0,0
    else: return ret[0],ret[1]


def plot_weighted_2d(x,y,weights,disc,bins,name):

    plt.plot(x,y,'.r')
    plt.xlabel('x')
    plt.ylabel('y')

    # Estimate the 2D histogram
    nbins = bins
    H1, xedges, yedges = np.histogram2d(x,y,bins=nbins)
    print 'The shap of h1 is ' + str( H1.shape)
    H, xedges, yedges = np.histogram2d(x,y,bins=nbins,weights = weights)#np.zeros(len(y)) + 1/float(len(y)),normed=True)
    H1[H==0] = 1
    H[H==0] = 1

    H = H/H1
    # H needs to be rotated and flipped

    H = np.rot90(H)
    H = np.flipud(H)

    # Mask zeros
    Hmasked = np.ma.masked_where(H==0,H) # Mask pixels with a value of zero

    # Plot 2D histogram using pcolor
    plt.pcolormesh(xedges,yedges,Hmasked)
    output = open('../output/pkl/eff.pkl', 'wb')

    pkl.dump({'eff_tuple':(xedges,yedges,Hmasked)},output)
    plt.xlabel('$\\Delta$ R(e-$\\gamma$)',fontsize = 20)
    plt.ylabel('$\\gamma$ Energy in MeV',fontsize = 20)
    plt.title('$\\epsilon$(disc > %s), NEntries = %s' %(disc,len(x)),fontsize=20)
    cbar = plt.colorbar()
    cbar.ax.set_ylabel('Efficiency for Bkg. Rejection')
    plt.savefig(name)
    plt.clf()


def plot_sig_bkg_dists(x_titles,t_vars,t_target):
    assert(len(x_titles) == len(t_vars[0])-2 )
    for line in range(2,len(t_vars[0])-1):
        bin_min = min( np.min(t_vars[:,line][t_target == 0]),np.min(t_vars[:,line][t_target == 1]) )
        bin_max = max( np.max(t_vars[:,line][t_target == 0]),np.max(t_vars[:,line][t_target == 1]) )
        weights_1 = np.ones_like(t_vars[:,line][t_target == 0])/len(t_vars[:,0][t_target == 0]);
        weights_2 = np.ones_like(t_vars[:,line][t_target == 1])/len(t_vars[:,0][t_target == 1]);
        bins    = np.arange(bin_min,bin_max,(bin_max-bin_min)/60.0)
        plt.hist(t_vars[:,line][t_target == 0],color = 'blue',alpha=.5, weights = weights_1,bins = bins,label = 'Background Events')
        plt.hist(t_vars[:,line][t_target == 1],color='red',alpha=.5,weights=weights_2,bins = bins,label = 'Signal Events')#bins=np.arange(np.min(v_vars[:,0][v_target == 1]),np.max(v_vars[:,0][v_target == 1]),30))
        plt.title('Geant4 Sim.',fontsize=16,loc='left')
        plt.ylabel('Normalized Dist.',fontsize = 20)
        plt.xlabel(x_titles[line-2],fontsize = 20)
        plt.legend(loc='upper right')
        plt.savefig('../output/plots/feature_%s.png' %(line))

        plt.clf()


def plot_binned_rocs(preds_1_total_bkg,v_target_total_bkg,preds_1_total_sig,v_target_total_sig,sep_total_bkg):
    colors = ['blue','red','yellow','green','cyan']
    labels = ['$\\Delta R$  $\epsilon $ [0,2.5] mm','$\\Delta R$  $\epsilon $ [2.5,5] mm','$\\Delta R$  $\epsilon $ [5.,7.5] mm',
              '$\\Delta R$  $\epsilon $ [7.5,10] mm']

    plt.title('Receiver Operating Characteristic',fontsize=24)
    false_positive_rate, true_positive_rate, thresholds = roc_curve(np.append(v_target_total_bkg,v_target_total_sig),np.append(preds_1_total_bkg,preds_1_total_sig))
    plt.plot(false_positive_rate, true_positive_rate, color = colors[-1],
             label='Inclusive AUC = %0.10E  ' % ((Decimal( auc(false_positive_rate, true_positive_rate)))))


    for id,sep in enumerate([0,2.5,5.,7.5]):
            inds        = [id2 for id2,x in enumerate(sep_total_bkg) if sep < x < sep+2.5]
            preds_bin   = np.append(preds_1_total_bkg[inds],preds_1_total_sig)
            target_bin  = np.append(v_target_total_bkg[inds],v_target_total_sig)
            plt.title('Receiver Operating Characteristic',fontsize=24)
            false_positive_rate, true_positive_rate, thresholds = roc_curve(target_bin,preds_bin)
            plt.plot(false_positive_rate, true_positive_rate, color = colors[id],
                     label='AUC (%s)  = %0.10E  ' % (labels[id],(Decimal( auc(false_positive_rate, true_positive_rate)))))

    plt.legend(loc='lower right')
    plt.savefig('../output/plots/ROC.png')
    plt.clf()


def sweep(preds_1_total_sig,mis_tag = .05):
    for cut in np.array(range(10000)).astype(float)/10000:
        mistag = float(len(preds_1_total_sig[preds_1_total_sig < cut]))/len(preds_1_total_sig)
        if mistag > mis_tag: break
    return cut

