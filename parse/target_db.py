import numpy as np
from pymongo import MongoClient

reset_db    = True
client      = MongoClient()
db          = client['hgcal']

if reset_db:
    db.target_magnetic_field.remove()
for file in ['brem_t0','brem_t1','brem_t2','brem_t3','brem_t4','brem_t5','brem_t6','brem_t7']:
    in_name  = "../data/target_files/%s.txt" %(file)
    lines = [line.rstrip('\n') for line in open(in_name)]
    events,temps = [],[]
    prev_run = 0
    for line in lines[3:len(lines)-2]:
        split =  line.split("*")

        if prev_run != split[-2]:
            if len(temps) > 0:
                db.target_magnetic_field.insert({'particles': temps})
                temps = []
        if float(split[2]) > 0:
            temps.append({'pdg':11,'eng':float(split[2]),'theta':float(split[4]),'phi':float(split[6]),'run':split[8],'x_pos':float(split[9]),'y_pos':float(split[10]),'z_pos':float(split[11])})
        else:
            #temps.append([11,float(split[3]),float(split[5]),float(split[7]),float(split[-2])])
            temps.append({'pdg':22,'eng':float(split[3]),'theta':float(split[5]),'phi':float(split[7]),'run':split[8],'x_pos':float(split[9]),'y_pos':float(split[10]),'z_pos':float(split[11])})
        prev_run = split[-2]
    print np.array(events)
