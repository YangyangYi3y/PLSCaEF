
clear
% Component correlation Significance
ProjectFolder = '/GPFS/cuizaixu_lab_permanent/yiyangyang/projects/ABCD/abcd_ef/PLSca';
PLSca_Folder = [ProjectFolder,'/results/ef_task23/combat/2fold'];

%actual corr
Res_Cell_dir = dir([PLSca_Folder,'/RandomCV_101Repeats_RegressCovariates_All_2Fold/Time*/Res_NFold.mat']);
for i = 1:length(Res_Cell_dir)
    Res_Cell{i} = [Res_Cell_dir(i).folder,'/',Res_Cell_dir(i).name];
end
for i = 1:length(Res_Cell)
    tmp = load(Res_Cell{i});
    Corr_Actual(i, :) = tmp.Mean_Corr;
end
media_Corr_Actual=  median(Corr_Actual)';

% Permutation corr
Permutation_Cell_dir = dir([PLSca_Folder,'/RandomCV_RegressCovariates_All_2Fold_Permutation/Time*/Res_NFold.mat']);
for i = 1:length(Permutation_Cell_dir)
    Permutation_Cell{i} = [Permutation_Cell_dir(i).folder,'/',Permutation_Cell_dir(i).name];
end
for i = 1:length(Permutation_Cell)
  tmp = load(Permutation_Cell{i});
  Corr_Permutation(i, :) = tmp.Mean_Corr;
end
for i = 1:size(Corr_Actual,2)
    Corr_Sig(i) = length(find(Corr_Permutation(:,i) > media_Corr_Actual(i))) / length(Permutation_Cell);
end

% Average covariance
Cov_cell_dir = dir([PLSca_Folder,'/RandomCV_101Repeats_RegressCovariates_All_2Fold/Time*/Fold_*'])
for i = 1:length(Res_Cell_dir)
    Cov_cell{i} = [Cov_cell_dir(i).folder,'/',Cov_cell_dir(i).name];
end
for i = 1:length(Cov_cell)
  tmp = load(Cov_cell{i});
  CovarianceExplained(i, :) = tmp.Covariances;
end
CovarianceExplained_Median = median(CovarianceExplained);
save([PLSca_Folder, '/Comp_corr','.mat'], 'Corr_Actual', 'Corr_Permutation', 'media_Corr_Actual', 'Corr_Sig','CovarianceExplained', 'CovarianceExplained_Median');

