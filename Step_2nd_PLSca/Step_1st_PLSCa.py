import scipy.io as sio
import numpy as np
import os
import sys
sys.path.append('/code/Functions');
import PLSca_CZ_Random_RegressCovariates

ProjectFolder = '';
EF_Mat=sio.loadmat(ProjectFolder + '/EF_data.mat');
EF_Data = EF_Mat['EF_Data'];
Conn_Mat = sio.loadmat(ProjectFolder + '/conn_data_reorder.mat');
Conn_Data = Conn_Mat['Conn_Data_ReOrder'];
Covariates_Mat = sio.loadmat(ProjectFolder + '/Covariates.mat');
subjectnum = len(EF_Data);
Covariates = np.zeros((subjectnum,3));
Covariates[:,0] = Covariates_Mat['Sex'];
Covariates[:,1] = Covariates_Mat['Age'];
Covariates[:,2] = Covariates_Mat['FD']

Components_Number = 23
Fold_Quantity = 2
CVRepeatTimes = 101
Permutation_time = 1000


ResultantFolder = ProjectFolder +'/RandomCV_101Repeats_RegressCovariates_All_' + str(Fold_Quantity)+'Fold'
PLSca_CZ_Random_RegressCovariates.PLSca_KFold_RandomCV_MultiTimes(Conn_Data, EF_Data, Covariates, Components_Number, Fold_Quantity, CVRepeatTimes, ResultantFolder, 0)

ResultantFolder = ProjectFolder +'/RandomCV_RegressCovariates_All_' + str(Fold_Quantity)+'Fold_Permutation'
PLSca_CZ_Random_RegressCovariates.PLSca_KFold_RandomCV_MultiTimes(Conn_Data, EF_Data, Covariates, Components_Number, Fold_Quantity, Permutation_time , ResultantFolder, 1)


