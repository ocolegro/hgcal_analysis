#!/bin/bash
source /afs/cern.ch/user/o/ocolegro/standalone_sim/PFCalEE/g4env.sh
cp /afs/cern.ch/user/o/ocolegro/standalone_sim/PFCalEE/git_hexaV02-01-01/version_30/model_0/e-/BOFF//et_30/eta_10000.000//run_19//g4steer.mac .
PFCalEE g4steer.mac 30 0 10000.000000 1.75,1.75,1.75,1.75,1.75,2.8,2.8,2.8,2.8,2.8,4.2,4.2,4.2,4.2,4.2 1,1,1,1,1,2.1,2.1,2.1,2.1,2.1,4.4,4.4,4.4,4.4  | tee g4.log
mv PFcal.root HGcal__version30_model0_BOFF_et30_eta10000.000_run19.root
localdir=`pwd`
echo "--Local directory is " $localdir >> g4.log
ls * >> g4.log
grep "alias eos=" /afs/cern.ch/project/eos/installation/cms/etc/setup.sh | sed "s/alias /export my/" > eosenv.sh
source eosenv.sh
$myeos mkdir -p /afs/cern.ch/work/o/ocolegro/test/githexaV02-01-01/e-
cmsStage -f HGcal__version30_model0_BOFF_et30_eta10000.000_run19.root /afs/cern.ch/work/o/ocolegro/test/githexaV02-01-01/e-/HGcal__version30_model0_BOFF_et30_eta10000.000_run19.root
if (( "$?" != "0" )); then
echo " --- Problem with copy of file PFcal.root to EOS. Keeping locally." >> g4.log
else
eossize=`$myeos ls -l /afs/cern.ch/work/o/ocolegro/test/githexaV02-01-01/e-/HGcal__version30_model0_BOFF_et30_eta10000.000_run19.root | awk '{print $5}'`
localsize=`ls -l HGcal__version30_model0_BOFF_et30_eta10000.000_run19.root | awk '{print $5}'`
if (( "$eossize" != "$localsize" )); then
echo " --- Copy of sim file to eos failed. Localsize = $localsize, eossize = $eossize. Keeping locally..." >> g4.log
else
echo " --- Size check done: Localsize = $localsize, eossize = $eossize" >> g4.log
echo " --- File PFcal.root successfully copied to EOS: /afs/cern.ch/work/o/ocolegro/test/githexaV02-01-01/e-/HGcal__version30_model0_BOFF_et30_eta10000.000_run19.root" >> g4.log
rm HGcal__version30_model0_BOFF_et30_eta10000.000_run19.root
fi
fi
echo "--deleting core files: too heavy!!"
rm core.*
cp * /afs/cern.ch/user/o/ocolegro/standalone_sim/PFCalEE/git_hexaV02-01-01/version_30/model_0/e-/BOFF//et_30/eta_10000.000//run_19//
echo "All done"
