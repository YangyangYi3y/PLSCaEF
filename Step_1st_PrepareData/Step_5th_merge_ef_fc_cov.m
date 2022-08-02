Conn_Mat = load('D:\projects\ABCD\func_baseline_10min\fc_combatdata.mat');
Conn_Data = Conn_Mat.combat_data;
Conn_Data = Conn_Data';
Conn_SubjectKey = Conn_Mat.ID{1};
ef = importdata('D:\projects\ABCD\EF\data\ef_task23.csv');
EF_Data = ef.data;
EF_ItemName = ef.textdata(1,3:end);
EF_SubjectKey = ef.textdata(2:end, 2);
Cov_data= importdata('D:\projects\ABCD\func_baseline_10min\covariables.txt');
Age_All = Cov_data.data(:, 1);
Sex_All = Cov_data.data(:, 2);
FD_All = Cov_data.data(:, 3);
Cov_SubjectKey = Cov_data.textdata(2:end, 1);

for i = 1:length(EF_SubjectKey)
  ind = find(strcmp(Conn_SubjectKey,EF_SubjectKey{i}));
  Conn_Data_ReOrder(i, :) = Conn_Data(ind, :);
  ind = find(strcmp(Cov_SubjectKey, EF_SubjectKey{i}));
  Age(i) = Age_All(ind(1));
  Sex(i) = Sex_All(ind(1));
  FD(i) = FD_All(ind(1));
end

save([ProjectFolder '/conn_data_reorder.mat'], 'Conn_Data_ReOrder');
save([ProjectFolder '/EF_data.mat'], 'EF_Data', 'EF_ItemName', 'EF_SubjectKey');
save([ProjectFolder '/Covariates.mat'], 'Age', 'Sex','FD');
