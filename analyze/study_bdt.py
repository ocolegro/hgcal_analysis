import sys
sys.path.insert(0,"../")
from sklearn.preprocessing import StandardScaler
from sklearn.cross_validation import StratifiedKFold
import xgboost as xgb
import numpy as np
import utils as u
from pymongo import MongoClient
from copy import copy


n_merged,n_unmerged,normed,max_sep,min_eng,max_eng  = 500000,500000,0,20,2000,4000
##Load up specified number of events from MongoDB
client      = MongoClient()
db          = client['hgcal']

if n_merged == -1:
    db_merged   = db.ele_gamma_finer.find(no_cursor_timeout=True).batch_size(10)
else:
    db_merged   = db.ele_gamma_finer.find(no_cursor_timeout=True).limit(n_merged)

if n_unmerged == -1:
    db_unmerged = db.electrons_finer.find(no_cursor_timeout=True).batch_size(10)
else:
    db_unmerged = db.electrons_finer.find(no_cursor_timeout=True).limit(n_unmerged)

merged_evts,unmerged_evts =  u.prep_bdt(u.load_event_array(db_merged),True), u.prep_bdt(u.load_event_array(db_unmerged),False)


'''


n_merged,n_unmerged,normed,max_sep,min_eng,max_eng  = 10000,10000,0,20,2000,4000
##Load up specified number of events from MongoDB
client      = MongoClient()
db          = client['hgcal']

if n_merged == -1:
    db_merged   = db.ele_gamma.find(no_cursor_timeout=True).batch_size(10)
else:
    db_merged   = db.ele_gamma.find(no_cursor_timeout=True).limit(n_merged)

if n_unmerged == -1:
    db_unmerged = db.electrons.find(no_cursor_timeout=True).batch_size(10)
else:
    db_unmerged = db.electrons.find(no_cursor_timeout=True).limit(n_unmerged)

merged_evts,unmerged_evts = u.prep_bdt(u.load_event_array(db_merged)),u.prep_bdt(u.load_event_array(db_unmerged),True)

for id,event in enumerate(db_merged):
    print (event)


print len(merged_evts),len(unmerged_evts),u.load_event_array(db_merged)


'''

##Build the important quantities up.
separation                          = copy(merged_evts[:,0])
unmerged_evts[:,0],merged_evts[:,0] = np.zeros(len(unmerged_evts)) + 1,np.zeros(len(merged_evts))
data,separation                     = np.vstack((unmerged_evts,merged_evts)),np.append(unmerged_evts[:,0],separation)

##Scaling helps the learner
if normed == True:
    data[:,range(2,len(unmerged_evts[0]))] = StandardScaler().fit_transform(copy(data[:,range(2,len(unmerged_evts[0]))]))


##Begin Stratified CV, I suggest you look into Sklearn's documentation to understand how this yields more statistics~
folds = StratifiedKFold(data[:,0], n_folds=2, shuffle=True, random_state=1)
counter = 0
sep_total_bkg,sep_total_sig,eng_total_bkg,eng_total_sig,weights_total,v_target_total_bkg,v_target_total_sig,preds_1_total_bkg,preds_1_total_sig = [],[],[],[],[],[],[],[],[]

for data_t,data_v in folds:
    t_vars,  v_vars    = data[data_t][:,range(2,len(data[0]))],data[data_v][:,range(2,len(data[0]))]
    t_target,v_target  = data[data_t][:,0], data[data_v][:,0]
    d_train,d_test     = xgb.DMatrix(t_vars,  label = t_target),xgb.DMatrix(v_vars,  label = v_target)


    param_1,num_round  = {'max_depth':3,'eta':.1, 'silent':1,'objective':'binary:logistic','eval_metric':'auc','nthread':6},300
    bst_1              = xgb.train(param_1,d_train,num_round)
    preds_1            = bst_1.predict(d_test)


    sep_total_bkg      = np.append(sep_total_bkg,separation[data_v][v_target==0])
    sep_total_sig      = np.append(sep_total_sig,separation[data_v][v_target==1])

    eng_total_bkg      = np.append(eng_total_bkg,data[data_v][:,1][v_target==0])
    eng_total_sig      = np.append(eng_total_sig,data[data_v][:,1][v_target==1])


    v_target_total_bkg = np.append(v_target_total_bkg,v_target[v_target==0])
    preds_1_total_bkg  = np.append(preds_1_total_bkg,preds_1[v_target==0])
    v_target_total_sig = np.append(v_target_total_sig,v_target[v_target==1])
    preds_1_total_sig  = np.append(preds_1_total_sig,preds_1[v_target==1])

cut_        = u.sweep(preds_1_total_sig)
weights_bkg = np.array(preds_1_total_bkg<cut_).astype(int)

weights_bkg     = weights_bkg[[id for id,x in enumerate(eng_total_bkg) if min_eng < x < max_eng]]
sep_total_bkg   = sep_total_bkg[[id for id,x in enumerate(eng_total_bkg) if min_eng < x < max_eng]]
eng_total_bkg   = np.array([x for x in eng_total_bkg if min_eng < x < max_eng])

weights_bkg     = weights_bkg[[id for id,x in enumerate(sep_total_bkg) if  x < max_sep]]
eng_total_bkg   = eng_total_bkg[[id for id,x in enumerate(sep_total_bkg) if x < max_sep]]
sep_total_bkg   = np.array([x for x in sep_total_bkg if x < max_sep])


u.plot_weighted_2d(sep_total_bkg,eng_total_bkg,weights_bkg,cut_,50,'../output/plots/eff2')
u.plot_weighted_2d(sep_total_bkg,eng_total_bkg,weights_bkg,cut_,10,'../output/plots/eff')

x_titles = ['Number of Hits'                      ,'Energy Deposited',  'Max Layer Depth',
            '68 % Containment Radius (Layers 3-7)','90 %Containment Radius (Layers 3-7)',
            '68 % Containment Radius (Layers 3-27)','90 % Containment Radius (Layers 3-27)',
            'Fraction of Hits (Layer 3)','Fraction of Hits (Layer 4)','Fraction of Hits (Layer 5)',
            'Fraction of Energy (Layer 3)','Fraction of Energy (Layer 4)','Fraction of Energy (Layer 5)']

u.plot_sig_bkg_dists(x_titles,t_vars,t_target)

u.plot_binned_rocs(preds_1_total_bkg,v_target_total_bkg,preds_1_total_sig,v_target_total_sig,sep_total_bkg)


