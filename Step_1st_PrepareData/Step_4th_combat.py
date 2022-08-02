
from neuroCombat import neuroCombat
import pandas as pd
import mat73
import scipy.io as sio

data = mat73.loadmat(r'D:\projects\ABCD\func_baseline_10min\fcm.mat')
fc = data['fc'].T

covars = pd.read_csv(r'D:\projects\ABCD\func_baseline_10min\cov.txt')
cotegorical_cols = ['sex']
continuous_col = ['age']

# To specify the name of the variable that encodes for the scanner/batch covariate:
batch_col = 'batch'
#Harmonization step:
fc_combat = neuroCombat(dat=fc,
                        covars = covars,
                        batch_col=batch_col,
                        categorical_cols = cotegorical_cols,
                        continuous_cols=continuous_col,
                        mean_only = True)
fc_combatdata = {'combat_data':fc_combat['data']}
sio.savemat(r'D:\projects\ABCD\func_baseline_10min\fc_combatdata.mat', fc_combatdata)


