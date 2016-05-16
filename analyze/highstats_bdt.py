import sys
sys.path.insert(0,"../")
from sklearn.preprocessing import StandardScaler,Imputer
from sklearn.cross_validation import StratifiedKFold
import xgboost as xgb
import numpy as np
import utils as u
from pymongo import MongoClient
from copy import copy


n_train,n_test,normed,max_sep,min_eng,max_eng  = 50000,1000000,0,20,2000,4000
##Load up specified number of events from MongoDB
client          = MongoClient()
db              = client['hgcal']

merged_evt      = db.ele_gamma_lowgran.find(no_cursor_timeout=True).limit(n_train)
unmerged_evt    = db.electrons_lowgran.find(no_cursor_timeout=True).limit(n_train)
merged_test     = db.ele_gamma_lowgran.find(no_cursor_timeout=True).skip(n_train).limit(n_test)


merged_evts,unmerged_evts =  u.prep_bdt(u.load_event_array(merged_evt) ,True,False), u.prep_bdt(u.load_event_array(unmerged_evt) ,False,False)
test                      =  u.prep_bdt(u.load_event_array(merged_test),True,False)
print len(test)
    
##Build the important quantities up.
sep_evts,unmerged_evts[:,0],merged_evts[:,0] \
                                    = copy(merged_evts[:,0]), np.zeros(len(unmerged_evts)) + 1, np.zeros(len(merged_evts))
evts,sep_evts                       = np.vstack((unmerged_evts,merged_evts)), np.append(unmerged_evts[:,0],sep_evts)

sep_test                            = copy(test[:,0])
test[:,0]                           = np.zeros(len(test)) + 1


##Scaling helps the learner
if normed == True:
    evts[:,range(2,len(unmerged_evts[0]))] = StandardScaler().fit_transform(copy(evts[:,range(2,len(unmerged_evts[0]))]))
    test[:,range(2,len(test[0]))]          = StandardScaler().fit_transform(copy(test[:,range(2,len(test[0]))]))

else:
    evts[:,range(2,len(unmerged_evts[0]))] = copy(evts[:,range(2,len(unmerged_evts[0]))])
    test[:,range(2,len(test[0]))]          = copy(test[:,range(2,len(test[0]))])

##Begin Stratified CV, I suggest you look into Sklearn's documentation to understand how this yields more statistics~
inds                        = np.random.permutation((len(evts)))
fold_1,fold_2               = inds[0:int(.8*len(evts))],inds[int(.8*len(evts)):]

t_vars,  t_target, v_vars, v_target \
                            = evts[:,range(2,len(evts[0]))][fold_1],evts[:,0][fold_1],evts[:,range(2,len(evts[0]))][fold_2],evts[:,0][fold_2]

tt_vars, tt_target          = test[:,range(2,len(test[0]))],test[:,0]

d_train, d_vars, d_test    \
                            = xgb.DMatrix(t_vars,  label = t_target), xgb.DMatrix(v_vars,  label = v_target), xgb.DMatrix(tt_vars,  label = tt_target)

param_1,num_round           = {'max_depth':3,'eta':.1, 'silent':1,'objective':'binary:logistic','eval_metric':'auc','nthread':6},50
bst_1                       = xgb.train(param_1,d_train,num_round)

preds_valid                 = bst_1.predict(d_vars)
preds_test                  = bst_1.predict(d_test)
print len(preds_valid),len(preds_valid[v_target==1])
cut_                        = u.sweep(preds_valid[v_target==1])



eng                = test[:,1]
weights_bkg        = np.array(preds_test < cut_).astype(int)


weights_bkg        = weights_bkg[[id for id,x in enumerate(eng) if min_eng < x < max_eng]]
sep_total_bkg      = sep_test[[id for id,x in enumerate(eng) if min_eng < x < max_eng]]
eng_total_bkg      = np.array([x for x in eng if min_eng < x < max_eng])


weights_bkg        = weights_bkg[[id for id,x in enumerate(sep_total_bkg) if  x < max_sep]]
eng_total_bkg      = eng_total_bkg[[id for id,x in enumerate(sep_total_bkg) if x < max_sep]]
sep_total_bkg      = np.array([x for x in sep_total_bkg if x < max_sep])


u.plot_weighted_2d(sep_total_bkg,eng_total_bkg,weights_bkg,cut_,50,'../output/plots/eff2_full_low_gran')
u.plot_weighted_2d(sep_total_bkg,eng_total_bkg,weights_bkg,cut_,10,'../output/plots/eff_full_low_gran')