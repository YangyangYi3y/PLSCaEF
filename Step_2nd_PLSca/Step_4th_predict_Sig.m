clear
%using brain component to predict 2 years latter behaviral component
fc = load('D:\projects\ABCD\EF\PLS\results\ef_task23\conn_data_reorder.mat');
fc_data = fc.Conn_Data_ReOrder;
ef= load('D:\projects\ABCD\EF\PLS\results\ef_task23\EF_data.mat');
subjectkey =ef.EF_SubjectKey;
ef_data = ef.EF_Data;

component = 2;
ProjectFolder = 'D:\projects\ABCD\EF\PLS\results\ef_task23\combat\2fold\';
PLSca_Folder = load([ProjectFolder,'Weight_sig_Comp',num2str(component),'.mat']);
Brain_Weight_actual = PLSca_Folder.Brain_Weight_actual;
Brain_Weight_median = median(Brain_Weight_actual);
Behavior_weight_median = PLSca_Folder.Behavior_weight_median;
Behavior_results = ef_data*Behavior_weight_median';
Brain_result = fc_data*Brain_Weight_median';

%sort component score
[pre_data_s,pre_index] = sort(Behavior_results);
sk_preL20 = subjectkey(pre_index(1:round(length(pre_index)*0.5)));
sk_preH20=subjectkey(pre_index(round(length(pre_index)*0.5)+1:length(pre_index)));

%the cbcl beseline score difference between high and low 20% component score
cbcl_dsm5Baseline = importdata('D:\projects\ABCD\EF\data\cbcl_dsm5_baseline.csv');
dsm_name = cbcl_dsm5Baseline.textdata(1,2:end);
cb_b_data = cbcl_dsm5Baseline.data;
cb_b_subjectkey = cbcl_dsm5Baseline.textdata(2:end,1);
nan_index = find(sum(isnan(cb_b_data),2)~=0);
cb_b_subjectkey(nan_index)= [];
cb_b_data(nan_index,:)=[];
cb_b_subjectkey = erase(cb_b_subjectkey,'_');
cd_b_preL_subjectkey = intersect(sk_preL20,cb_b_subjectkey);
cd_b_beh_subjectkey = intersect(sk_preH20,cb_b_subjectkey);
for i = 1:size(cb_b_data,2)
    hist(cb_b_data(:,i),50)
    title(char(dsm_name(i)));
    filename = ['D:\projects\ABCD\EF\PLS\results\ef_task23\combat\2fold\',char(dsm_name(i)),'_orig_hist.png']
    saveas(gcf,filename)
    close gcf
end


%difference in factor score of psychiatry between high and low 20% component score
factor = importdata('D:\projects\ABCD\EF\data\factorscore_wx_1y.csv');
factor_data = factor.data;
factor_name = factor.textdata(1,2:end);
factor_sk = factor.textdata(2:end,1);
factor_sk = erase(factor_sk,'_');
nan_index = find(sum(isnan(factor_data),2)~=0);
factor_sk(nan_index)= [];
factor_data(nan_index,:)=[];
factor_preL_subjectkey = intersect(sk_preL20,factor_sk);
factor_preH_subjectkey = intersect(sk_preH20,factor_sk);

%histogram of factor score of psychiatry
for i = 1:size(factor_data,2)
    name = char(factor_name(i));
    h1 = histogram(factor_data(:,i));
    nbins = 50
    title(name)
    saveas(gcf,['D:\projects\ABCD\EF\PLS\results\ef_task23\combat\2fold\',name,'_1y.png'])
end

%reorder factor score and component score
f_ef_subjectkey = intersect(subjectkey,factor_sk);
for i = 1:length(f_ef_subjectkey)
  ind = find(strcmp(factor_sk,(f_ef_subjectkey{i})));
  factor_data_reorder(i, :) = factor_data(ind, :);
  ind =  find(strcmp(subjectkey,(f_ef_subjectkey{i})));
  Behavior_results_reorder(i, :) = Behavior_results(ind, :);
end
cor_ef_f = corr(Behavior_results_reorder,factor_data_reorder);

% factor score of high and low 20% participants of component score
for i = 1:length(factor_preL_subjectkey)
  ind = find(strcmp(factor_sk,(factor_preL_subjectkey{i})));
  factorL20_reorder(i, :) = factor_data(ind, :);
end
for i = 1:length(factor_preH_subjectkey)
  ind =  find(strcmp(factor_sk,(factor_preH_subjectkey{i})));
  factorH20_reorder(i, :) = factor_data(ind, :);
end
total_num =  length(factorL20_reorder)+length(factorH20_reorder);
factor20= zeros(total_num,5);
factor20(1:length(factorL20_reorder),1:4) = factorL20_reorder;
factor20(1:length(factorL20_reorder),5) = 1;
factor20(length(factorL20_reorder)+1:total_num,1:4) = factorH20_reorder;
factor20(length(factorL20_reorder)+1:total_num,5) = 0;

for i = 1:size(factorL20_reorder,2)
   [h(i),p(i),ci(i,:),stats(i)] = ttest2(factor20(:,i),factor20(:,5));
   h1 = histogram(factorL20_reorder);
   hold on
   h2 = histogram(factorH20_reorder);
   title(char(factor_name(i)))
   filename = ['D:\projects\ABCD\EF\PLS\results\ef_task23\combat\2fold\comp',num2str(component),'\',char(factor_name(i)),'_his_1y.png'];
   saveas(gcf,filename)
   close gcf
   boxplot(factor20(:,i),factor20(:,5))
   title(char(factor_name(i)))
   filename = ['D:\projects\ABCD\EF\PLS\results\ef_task23\combat\2fold\comp',num2str(component),'\',char(factor_name(i)),'_box_1y.png'];
   saveas(gcf,filename)
   close gcf
end

% the difference of component score between high and low 10% factor pychiatry scores
for i = 1:size(factor_data,2)
    name = char(factor_name(i));
    data = factor_data(:,i);
    [data1,index] = sort(data);
    sk_L10 = factor_sk(index(1:round(length(factor_sk)*0.5)));
    sk_H10 = factor_sk(index(round(length(factor_sk))*0.5+1:length(factor_sk)));
    skL10_factor_ef = intersect(subjectkey,sk_L10);
    %skh20_de_ef = intersect(subjectkey,depresssh20_sk);
    skH10_factor_ef = intersect(subjectkey,sk_H10);
    for j = 1:length(skL10_factor_ef)
        ind = find(strcmp(subjectkey,(skL10_factor_ef{j})));
        factorL10(j) = Behavior_results(ind);
    end
    for j = 1:length( skH10_factor_ef)
        ind = find(strcmp(subjectkey,(skH10_factor_ef{j})));
        factorH10(j) = Behavior_results(ind);
    end
    h1 = histogram(factorL10);
    nbins = 50
    hold on
    h2 = histogram(factorH10);
    nbins = 50
    title(name)
    %saveas(gcf,['D:\projects\ABCD\EF\PLS\results\ef_task23\combat\2fold\comp',num2str(component),'\',name,'10',char(factor_name(i)),'_1y.png'])
    close gcf
    total_sub = length(factorL10)+ length(factorH10);
    factor_score(1:length(factorL10),1) = factorL10;
    factor_score(1:length(factorL10),2) = 0;
    factor_score(length(factorL10)+1:total_sub,1) = factorH10;
    factor_score(length(factorL10)+1:total_sub,2) = 1;
    d = struct('a',factor_score(:,1),'b',factor_score(:,2))
    [h10(i),p10_2(i),ci10(i,:),stats10(i)]= ttest2(factorL10,factorH10);
    boxplot(d.a,d.b)
    title(['component',num2str(component), 'score between 10% high',name, 'and 10% low'])
    xlabel(name)
    ylabel(['component',num2str(component),'score'])
    %saveas(gcf,['D:\projects\ABCD\EF\PLS\results\ef_task23\combat\2fold\comp',num2str(component),'\',name,'10',char(factor_name(i)),'_box_1y.png'])
    close gcf
end

for i = 1:length(cd_b_preL_subjectkey)
  ind = find(strcmp(cb_b_subjectkey,(cd_b_preL_subjectkey{i})));
  cd_bl20_reorder(i, :) = cb_b_data(ind, :);
  ind =  find(strcmp(cb_b_subjectkey,(cd_b_beh_subjectkey{i})));
  cd_bh20_reorder(i, :) = cb_b_data(ind, :);
end
total_num =  length(cd_bl20_reorder)+length(cd_bh20_reorder);
cbcl_base = zeros(total_num,7);
cbcl_base(1:length(cd_bl20_reorder),1:6) = cd_bl20_reorder;
cbcl_base(1:length(cd_bl20_reorder),7) = 1;
cbcl_base(length(cd_bl20_reorder)+1:total_num,1:6) = cd_bh20_reorder;
cbcl_base(length(cd_bl20_reorder)+1:total_num,7) = 0;



for i = 1:size(cd_bl20_reorder,2)
   [h(i),p(i),ci(i,:),stats(i)] = ttest2(cbcl_base(:,i),cbcl_base(:,7),'Vartype','unequal');
   h1 = histogram(cd_bl20_reorder);
   hold on
   h2 = histogram(cd_bh20_reorder);
   title(char(dsm_name(i)))
   filename = ['D:\projects\ABCD\EF\PLS\results\ef_task23\combat\2fold\comp',num2str(component),'\',char(dsm_name(i)),'_his.png'];
   saveas(gcf,filename)
   close gcf
   boxplot(cbcl_base(:,i),cbcl_base(:,7))
   title(char(dsm_name(i)))
   filename = ['D:\projects\ABCD\EF\PLS\results\ef_task23\combat\2fold\comp',num2str(component),'\',char(dsm_name(i)),'_box.png'];
   saveas(gcf,filename)
   close gcf
end

for i = 1:size(cb_b_data,2)
    name = char(dsm_name(i));
    data = cb_b_data(:,i);
    [data1,index] = sort(data);
    num= length(find(data1 ==50));
    sk_50 = cb_b_subjectkey(index(1:num));
    sk_H10 = cb_b_subjectkey(index(length(cb_b_subjectkey)*0.9+1:length(cb_b_subjectkey)));
    sk50_cb_ef = intersect(subjectkey,sk_50 );
    %skh20_de_ef = intersect(subjectkey,depresssh20_sk);
    skh10_de_ef = intersect(subjectkey,sk_H10);
    for j = 1:length(sk50_cb_ef)
        ind = find(strcmp(subjectkey,(sk50_cb_ef{j})));
        cb50(j) = Behavior_results(ind);
    end
    for j = 1:length(skh10_de_ef)
        ind = find(strcmp(subjectkey,(skh10_de_ef{j})));
        cbh10(j) = Behavior_results(ind);
    end
    h1 = histogram(cb50);
    nbins = 50
    hold on
    h2 = histogram(cbh10);
    nbins = 50
    title(name)
    saveas(gcf,['D:\projects\ABCD\EF\PLS\results\ef_task23\combat\2fold\comp',num2str(component),'\',name,'50_10_comp',num2str(component),'.png'])
    close gcf
    total_sub = length(cb50)+ length(cbh10);
    cb(1:length(cb50),1) = cb50;
    cb(1:length(cb50),2) = 0;
    cb(length(cb50)+1:total_sub,1) = cbh10;
    cb(length(cb50)+1:total_sub,2) = 1;
    d = struct('a',cb(:,1),'b',cb(:,2))
    [h50(i),p50(i),ci50(i,:),stats50(i)]= ttest2(cb50,cbh10);
    boxplot(d.a,d.b)
    title(['component',num2str(component), 'score between 10% high',name, 'and 50'])
    xlabel(name)
    ylabel(['component',num2str(component),'score'])
    saveas(gcf,['D:\projects\ABCD\EF\PLS\results\ef_task23\combat\2fold\comp',num2str(component),'\',name,'50_10_comp',num2str(component),'_box.png'])
    close gcf
end

behav= load('D:\projects\ABCD\EF\PLS\results\ef_task23\EF_data.mat');
fc_subjectkey = behav.EF_SubjectKey;
two_year_ef = importdata('D:\projects\ABCD\EF\data\efz_task23_2year.csv');
two_year_ef_data = two_year_ef.data;
two_year_subjectkey = two_year_ef.textdata(2:end,2);
behav_result = two_year_ef_data  *Behavior_weight_median';
subjectkey = intersect(two_year_subjectkey,fc_subjectkey);
for i = 1:length(subjectkey)
  ind = find(strcmp(fc_subjectkey,subjectkey{i}));
  brain_result_reorder(i, :) =brain_result(ind, :);
  ind = find(strcmp(two_year_subjectkey,subjectkey{i}));
  behav_result_reorder(i, :) = behav_result (ind, :);
end
[r,p] = corr(brain_result_reorder,behav_result_reorder)
save([ProjectFolder '/brain_2_2year_surveyComponent',num2str(component),'.mat'], 'brain_result_reorder','behav_result_reorder', 'r','p');

%using baseline bahavioral component to predict cbcl diagnose at baseline 
for i = 1:length(cd_b_ef_subjectkey)
  ind = find(strcmp(cb_b_subjectkey,(cd_b_ef_subjectkey{i})));
  cd_b_reorder(i, :) =cb_b_data(ind, :);
  ind =  find(strcmp(subjectkey,(cd_b_ef_subjectkey{i})));
  behavb_reorder(i, :) =Behavior_results(ind, :);
end
for i = 1:size(cd_b_reorder,2)
    [r,p] = corr(cd_b_reorder(:,i),behavb_reorder);
    r_b(i,1) = r;
    r_b(i,2) = p;
end

%using baseline brain component to predict cbcl diagnose at baseline 
for i = 1:length(cd_b_ef_subjectkey)
  ind = find(strcmp(cb_b_subjectkey,(cd_b_ef_subjectkey{i})));
  cd_b_reorder(i, :) =cb_b_data(ind, :);
  ind =  find(strcmp(subjectkey,(cd_b_ef_subjectkey{i})));
  brain_reorder(i, :) = Brain_result(ind, :);
end
for i = 1:size(cd_b_reorder,2)
    [r,p] = corr(cd_b_reorder(:,i),brain_reorder);
    r_b(i,1) = r;
    r_b(i,2) = p;
end

cbcl_dsm53year = importdata('D:\projects\ABCD\EF\data\cbcl_dsm5_3year.csv');
cd_3_data = cbcl_dsm53year.data;
nan_index = find(sum(isnan(cd_3_data),2)~=0);
cd_3_data(nan_index,:)=[];
cd_3_subjectkey = cbcl_dsm53year.textdata(2:end,1);
cd_3_subjectkey = erase(cd_3_subjectkey,'_');
cd_3_subjectkey(nan_index)=[];
cd_3_ef_subjectkey = intersect(cd_3_subjectkey,subjectkey);
for i = 1:length(cd_3_ef_subjectkey)
  ind = find(strcmp(cd_3_subjectkey,(cd_3_ef_subjectkey{i})));
  cd_3_reorder(i, :) =cd_3_data(ind, :);
  ind =  find(strcmp(subjectkey,(cd_3_ef_subjectkey{i})));
  behav3_reorder(i, :) =Behavior_results(ind, :);
end
for i = 1:size(cd_3_reorder,2)
    [r,p] = corr(cd_3_reorder(:,i),behav3_reorder);
    r_3(i,1) = r;
    r_3(i,2) = p;
end


fc = load('D:\projects\ABCD\EF\PLS\results\ef_task23\conn_data_reorder.mat');
fc_data = fc.Conn_Data_ReOrder;

fc_subjectkey = behav.EF_SubjectKey;
two_year_ef = importdata('D:\projects\ABCD\EF\data\efz_task23_2year.csv');
two_year_ef_data = two_year_ef.data;
cbcl_subjectkey = two_year_ef.textdata(2:end,2);
behav_result = two_year_ef_data *Behavior_weight_median';
subjectkey = intersect(two_year_subjectkey,fc_subjectkey);

[r,p] = corr(brain_result_reorder,behav_result_reorder)
save([ProjectFolder '/brain_2_2year_surveyComponent',num2str(component),'.mat'], 'brain_result_reorder','behav_result_reorder', 'r','p');

surveyBaseline = importdata("D:/projects/ABCD/EF/data/task23baseline.csv");
survey2years = importdata("D:/projects/ABCD/EF/data/task232years.csv");
surveyBdata = surveyBaseline.data;
surveyBsk = surveyBaseline.textdata(2:end,:);
survey2data = survey2years.data;
survey2sk = survey2years.textdata(2:end,:);
surveyBsk = erase(surveyBsk,',');
surveyBsk = erase(surveyBsk,'_');
survey2sk = erase(survey2sk,'_');
survey_sk = intersect(surveyBsk,survey2sk);
survey_sk(cellfun(@isempty,survey_sk))=[];
for i = 1:length(survey_sk)
  ind = find(strcmp(surveyBsk,(survey_sk{i})));
  surveyBdata_reorder(i, :) =surveyBdata(ind, :);
  ind =  find(strcmp(survey2sk,(survey_sk{i})));
  survey2data_reorder(i, :) =survey2data(ind, :);
end

sur_pre = corr(surveyBdata_reorder, survey2data_reorder);
heatmap(sur_pre)
B=[];
[m,n]=size(sur_pre);
for i=1:m
    B=[B;sur_pre(i,[1:i-1 i+1:n])];
end


ksad = importdata('D:\projects\ABCD\EF\data\ksads\ksads_diagnose_baseline_present2.txt');
health = importdata('D:\projects\ABCD\EF\data\ksads\ksads_health_subkey.csv');
ksad_diag = ksad.data(:,2:end);
ksad_name = ksad.textdata(1,4:end);
ksad_name = erase(ksad_name ,':');
ksad_sk = ksad.textdata(2:end,1);
ksad_sk  = erase(ksad_sk ,'_');
diag = zeros(4141,9);
for i =1:length(ksad_name)
    clear ksad_diagi_reorder
    clear behav_reorderi
    ksad_diagi = ksad_diag(:,i);
    nan_index = find(ismissing(ksad_diagi)==1);
    ksad_ski = ksad_sk;
    ksad_ski(nan_index)= [];
    ksad_diagi(nan_index,:)=[];
    inter_sk = intersect(subjectkey,ksad_ski);
    for j = 1:length(inter_sk)
        ind = find(strcmp(ksad_ski,(inter_sk{j})));
        ksad_diagi_reorder(j, :) = ksad_diagi(ind, :);
        ind =  find(strcmp(subjectkey,(inter_sk{j})));
        behav_reorderi(j, :) = Behavior_results(ind, :);
    end
    patient_indexi = find(ksad_diagi_reorder~=0);
    health_indexi = find(ksad_diagi_reorder==0);
    patient_numi = length(patient_indexi);
    patient_num1(i) = patient_numi;
    behav_health = behav_reorderi(health_indexi);
    behav_healthz = zscore(behav_health);
    index = find(behav_healthz>3);
    behav_healthz(index)=3;
    index = find(behav_healthz<-3);
    behav_healthz(index)= -3
    behav_patient = behav_reorderi(patient_indexi);
    behav_patientz = zscore(behav_patient);
    index = find(behav_patientz>3);
    behav_patientz(index)=3;
    index = find(behav_patientz<-3);
    behav_patientz(index)= -3
    %h1 = histogram(behav_health,'FaceColor','g');
    %nbins = 50
    %hold on
    %h2 = histogram(behav_patient,'FaceColor','r');
    %nbins = 20
    %h1.Normalization = 'probability';
    %h1.BinWidth = 0.25;
    %h2.Normalization = 'probability';
    %h2.BinWidth = 0.25;
    %name = [ksad_name{i},'patient',num2str(patient_num1)];
    %title(name)
    %saveas(gcf,['D:\projects\ABCD\EF\PLS\results\ef_task23\combat\2fold\comp',num2str(component),'\',name,'_brain.png'])
    %close gcf
    ksadi(1:patient_numi,1) = behav_patient;
    ksadi(1:patient_numi,2) = 1;
    ksadi(patient_numi+1:length(inter_sk),1) = behav_health;
    ksadi(patient_numi+1:length(inter_sk),2) = 0;
    d = struct('a',ksadi(:,1),'b',ksadi(:,2));
    %[Mean(i),SD(i),Z(i),p_value_one_tail(i),p_value_two_tails(i)] = mwwtest(behav_health',behav_patient');
    a = 'cohen'
    Effect(i) = computeCohen_d(behav_health,behav_patient);
    mwwi=mwwtest(behav_health',behav_patient'); 
    mwwip=mwwi.p;
    mwwp(i,1:2) = mwwip;
    [hc2n(i),pc2n(i),cic2n(i,:),statsc2n(i)]= ttest2(behav_health,behav_patient);
    diag(patient_indexi,i)=1;
    %boxplot(d.a,d.b)
    %title(['component',num2str(component)])
    %xlabel(name)
    %ylabel(['component',num2str(component),'score'])
    %saveas(gcf,['D:\projects\ABCD\EF\PLS\results\ef_task23\combat\2fold\comp',num2str(component),'\',name,'_brain_box.png'])
    %close gcf
end
    
    
factor_sk(nan_index)= [];
factor_data(nan_index,:)=[];
nan_index = find(sum(isnan(ksad_diag()),2)~=0);
factor_sk(nan_index)= [];
factor_data(nan_index,:)=[];
factor_preL_subjectkey = intersect(sk_preL20,factor_sk);
factor_preH_subjectkey = intersect(sk_preH20,factor_sk);


