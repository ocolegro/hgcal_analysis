# hgcal_analysis
Analysis repository for generating studies of the HGCAL
-- Requirements --
C++, Python, Root

Python Packages --
For Data
Pymongo: https://api.mongodb.com/python/current/api/pymongo/

For Analysis --
Matplotlib: http://matplotlib.org/
Scipy: https://www.scipy.org/
Numpy: http://www.numpy.org/

For Machine Learning --
XGBoost: https://github.com/dmlc/xgboost
SKLearn: http://scikit-learn.org/stable/

Workflow --

Step 1-- Use Geant4 and target_event_generator to create target ntuples.

Step 2-- Use Geant4 and PFCAL (https://github.com/pfs/PFCal) to generate HGCAL events.  (My personal modifications to this repository will be included here shortly).

Step 3-- Copy files locally

Step 4-- Convert files to MongoDB

Step 5-- Carry out analysis
