% Feature_Updata_and_sig
clear
component = 1;
PLSca_Folder = '/GPFS/cuizaixu_lab_permanent/yiyangyang/projects/ABCD/abcd_ef/PLSca/results/task23/combat/2fold';
Res_Cell_dir = dir([PLSca_Folder,'/RandomCV_101Repeats_RegressCovariates_All_2Fold/Time_*/Fold_*.mat']);
for i = 1:length(Res_Cell_dir)
    Res_Cell{i} = [Res_Cell_dir(i).folder,'/',Res_Cell_dir(i).name];
end
%compare the first CV with all other fold and CV, then correct
%the sign of weight to make all the weight having the same direction
tmp1 = load(Res_Cell{1});
Behavior_Weight1 = tmp1.Behavior_Weight(:, component);
for i = 1:length(Res_Cell)
  tmp = load(Res_Cell{i});
  corr = corrcoef(Behavior_Weight1,tmp.Behavior_Weight(:,component));
  corr = corr(1,2);
  if corr>0
    Behavior_Weight_actual(:,i) = tmp.Behavior_Weight(:,component);
    Brain_Weight_actual(:,i) = tmp.Brain_Weight(:,component);   
  else
    Behavior_Weight_actual(:,i) = -tmp.Behavior_Weight(:,component); 
    Brain_Weight_actual(:,i) = -tmp.Brain_Weight(:,component);    
  end
end

Brain_weight_median= median(Brain_Weight_actual');
Behavior_weight_median = median(Behavior_Weight_actual');


%Permutation
Res_Cell_dir = dir([PLSca_Folder,'/RandomCV_RegressCovariates_All_2Fold_Permutation/Time_*/Fold_*.mat']);
for i = 1:length(Res_Cell_dir)
    Res_Cell{i} = [Res_Cell_dir(i).folder,'/',Res_Cell_dir(i).name];
end
for i = 1:length(Res_Cell)
  tmp = load(Res_Cell{i});
  corr = corrcoef(Behavior_Weight1,tmp.Behavior_Weight(:,component));
  corr = corr(1,2);
  if corr>0
    Behavior_Weight_permu(:,i) = tmp.Behavior_Weight(:,component);
    Brain_Weight_permu(:,i) = tmp.Brain_Weight(:,component);     
  else
    Behavior_Weight_permu(:,i) = -tmp.Behavior_Weight(:,component);
    Brain_Weight_permu(:,i) = -tmp.Brain_Weight(:,component);   
  end
end

%calculate sig with permutation
for i = 1:length(Behavior_weight_median)
    Behavior_weight_sig(i) = length(find(abs(Behavior_weight_median(i))<abs(Behavior_Weight_permu(i,:))))./size(Behavior_Weight_permu,2);
end

% calculate the mapping matrix from fc vector to fc matrix
network = importdata('/GPFS/cuizaixu_lab_permanent/yiyangyang/projects/ABCD/func_baseline_10min/net_roi.csv');
roi_num = network.data;
network_name = network.textdata(2:end,1);
network_name_unique = unique(network_name);
fc_weight = zeros(352,352);
startindex = 1;
map = zeros(352,2);
for n =1:length(network_name_unique)
    network_namei = network_name_unique(n);
    index = find(strcmp(network_name,network_namei)==1);
    endindex = startindex + length(index)-1;
    roi_numi = roi_num(index);
    map(startindex:endindex,1) = roi_numi;
    map(startindex:endindex,2) = index;
    startindex = endindex+1;
end

% mapping from fc vector to fc matrix at roi level for actual brain weight 
Brain_Weight_actual = Brain_Weight_actual';
for s = 1:size(Brain_Weight_actual,1)    
    brain_weighti = Brain_Weight_actual(s,:);
    n = 1;
    fct = zeros(352,352);
    for j = 2:352
        for i = 1:352
            if j-1>=i
               fct(i,j)= brain_weighti(n);
               n = n+1;
            end
            i =i+1;   
        end
    end
    fct_1 = fct';
    corm= fct + fct_1; 
    for x = 1:352
        for y = 1:352
            new_i = map(x,1);
            new_j = map(y,1);
        fc_weight (x,y) =  corm(new_i,new_j) ;
        end
    end 
    fc_roi_weight_actual(s,:,:) = fc_weight;
end

% calculate the median roi weight across fold_number*repetition_number(202 for 2 fold)
 for i = 1: size(fc_roi_weight_actual,2)
     for j = 1:size(fc_roi_weight_actual,3)
        ij_t = fc_roi_weight_actual(:,i,j);
        fc_roi_MedianWeight (i,j) = median(ij_t);
     end
 end
 
% mapping from fc vector to fc matrix at roi level for permulation brain weight  
Brain_Weight_permu = Brain_Weight_permu';
 for s = 1:size(Brain_Weight_permu,1)    
    brain_weighti = Brain_Weight_permu(s,:);
    n = 1;
    fct = zeros(352,352);
    for j = 2:352
        for i = 1:352
            if j-1>=i
               fct(i,j)= brain_weighti(n);
               n = n+1;
            end
            i =i+1;   
        end
    end
    fct_1 = fct';
    corm= fct + fct_1; 
    for x = 1:352
        for y = 1:352
            new_i = map(x,1);
            new_j = map(y,1);
        fc_weight_per (x,y) =  corm(new_i,new_j) ;
        end
    end 
    fc_roi_weight_permu(s,:,:) = fc_weight_per;
end

 for i = 1: size(fc_roi_weight_actual,2)
     for j = 1:size(fc_roi_weight_actual,3)
         Brain_roi_weight_sig(i,j) = length(find(abs(fc_roi_MedianWeight(i,j))<abs(fc_roi_weight_permu(:,i,j))))./size(fc_roi_weight_permu,1);
     end
 end
 
% reduce roi level to networklevel
network_name0 = network.textdata(2:end,1);
network_name = unique(network_name0);
network_name(12)=[];
network_name(14) = cellstr('Subcortical');

%get the index for the roi in network 
index1 = [];
for i =1: length(network_name)
    network_namei = network_name(i);
    num= strmatch(network_namei,network_name0);
    index1 = [index1;num(end)];
end
index(1) = 1;
index(2:length(index1)+1)= index1;

%calculate network_level actual weight
fc_total_abs = zeros(14,14);
fc_total_p = zeros(14,14);
fc_total_n = zeros(14,14);
fc_net_weight_actual_sumn = zeros(202,14,14);
fc_net_weight_actual_abs =zeros(202,14,14);
fc_net_weight_actual_p =zeros(202,14,14) ;
fc_net_weight_actual_n =zeros(202,14,14); 
for s = 1:size(Brain_Weight_actual,1)  
    fc_weight = squeeze(fc_roi_weight_actual(s,:,:));
    for i = 1:length(index1)
        startindex_row = index(i);
        endindex_row = index(i+1);
        num_row = endindex_row-startindex_row+1;
        fc1 = fc_weight(startindex_row:endindex_row,:);
        for j = 1:length(index1)
            startindex_col = index(j);
            endindex_col= index(j+1);
            num_col = endindex_col - startindex_col+1;
            fc_neti = fc1(:,startindex_col:endindex_col);
            if i == j %within network, the diagonals are zeros, not count
                tnum = num_row  * (num_col-1);
            else % between network 
                tnum = num_row * (num_col);
            end
            fc_mean_sumn = sum(fc_neti,'ALL')/tnum;
            fc_mean_abs = sum(abs(fc_neti),'ALL')/tnum;
            index_p = find(fc_neti>0);
            fc_mean_p = sum(fc_neti(index_p),'ALL')/tnum;
            index_n = find(fc_neti<0);
            fc_mean_n = sum(fc_neti(index_n),'ALL')/tnum;
            fc_total_sumn(i,j) = fc_mean_sumn;
            fc_total_abs(i,j) = fc_mean_abs;
            fc_total_p(i,j) = fc_mean_p;
            fc_total_n(i,j) = fc_mean_n;
        end
        fc_net_weight_actual_sumn(s,:,:) = fc_total_sumn;
        fc_net_weight_actual_abs(s,:,:) = fc_total_abs;
        fc_net_weight_actual_p(s,:,:) = fc_total_p;
        fc_net_weight_actual_n(s,:,:) = fc_total_n;
    end
end


 for i = 1: size(fc_net_weight_actual_sumn,2)
     for j = 1:size(fc_net_weight_actual_sumn,3)
        ij_t = fc_net_weight_actual_sumn(:,i,j);
        fc_net_MedianWeight_sumn (i,j) = median(ij_t);
     end
 end
 for i = 1: size(fc_net_weight_actual_abs,2)
     for j = 1:size(fc_net_weight_actual_abs,3)
        ij_t = fc_net_weight_actual_abs(:,i,j);
        fc_net_MedianWeight_abs (i,j) = median(ij_t);
     end
 end
  for i = 1: size(fc_net_weight_actual_p,2)
     for j = 1:size(fc_net_weight_actual_p,3)
        ij_t = fc_net_weight_actual_p(:,i,j);
        fc_net_MedianWeight_p (i,j) = median(ij_t);
     end
  end
  for i = 1: size(fc_net_weight_actual_n,2)
     for j = 1:size(fc_net_weight_actual_n,3)
        ij_t = fc_net_weight_actual_n(:,i,j);
        fc_net_MedianWeight_n (i,j) = median(ij_t);
     end
 end
 

%calculate network_level permu weight
fc_total_per = zeros(14,14);
fc_total_p= zeros(14,14);
fc_total_n = zeros(14,14);
fc_net_weight_per_sumn =  zeros(2000,14,14);
fc_net_weight_per_abs =  zeros(2000,14,14);
fc_net_weight_per_p =  zeros(2000,14,14);
fc_net_weight_per_n =  zeros(2000,14,14);
%fc_net_weight_per = zeros(size(Brain_Weight_permu,1),length(index1),length(index1));
for s = 1:size(Brain_Weight_permu,1)
    fc_weight1 = squeeze(fc_roi_weight_permu(s,:,:));
    for i = 1:length(index1)
        startindex_row = index(i);
        endindex_row = index(i+1);
        num_row = endindex_row-startindex_row+1;
        fc1 = fc_weight1(startindex_row:endindex_row,:);
        for j = 1:length(index1)
            startindex_col = index(j);
            endindex_col= index(j+1);
            num_col = endindex_col - startindex_col+1;
            fc_neti = fc1(:,startindex_col:endindex_col);
            if i == j
                tnum = num_row  * (num_col-1);
            else
                tnum = num_row * (num_col);
            end
            fc_mean_sumn = sum(fc_neti,'ALL')/tnum;
            fc_mean_abs = sum(abs(fc_neti),'ALL')/tnum;
            index_p = find(fc_neti>0);
            fc_mean_p = sum(fc_neti(index_p),'ALL')/tnum;
            index_n = find(fc_neti<0);
            fc_mean_n = sum(fc_neti(index_n),'ALL')/tnum;
            fc_total_sumn(i,j) = fc_mean_sumn;
            fc_total_abs(i,j) = fc_mean_abs;
            fc_total_p(i,j) = fc_mean_p;
            fc_total_n(i,j) = fc_mean_n;
        end
        fc_net_weight_per_sumn(s,:,:) = fc_total_sumn;
        fc_net_weight_per_abs(s,:,:) = fc_total_abs;
        fc_net_weight_per_p(s,:,:) = fc_total_p;
        fc_net_weight_per_n(s,:,:) = fc_total_n;
    end
end

 for i = 1: size(fc_net_weight_actual_sumn,2)
     for j = 1:size(fc_net_weight_actual_sumn,3)
         Brain_net_weight_sig_sumn(i,j) = length(find(abs(fc_net_MedianWeight_sumn(i,j))<abs(fc_net_weight_per_sumn(:,i,j))))./size(fc_net_weight_per_sumn,1);
     end
 end
 for i = 1: size(fc_net_weight_actual_abs,2)
     for j = 1:size(fc_net_weight_actual_abs,3)
         Brain_net_weight_sig_abs(i,j) = length(find(abs(fc_net_MedianWeight_abs(i,j))<abs(fc_net_weight_per_abs(:,i,j))))./size(fc_net_weight_per_abs,1);
     end
 end
  for i = 1: size(fc_net_weight_actual_p,2)
     for j = 1:size(fc_net_weight_actual_p,3)
         Brain_net_weight_sig_p(i,j) = length(find(abs(fc_net_MedianWeight_p(i,j))<abs(fc_net_weight_per_p(:,i,j))))./size(fc_net_weight_per_p,1);
     end
  end
  for i = 1: size(fc_net_weight_actual_n,2)
     for j = 1:size(fc_net_weight_actual_n,3)
         Brain_net_weight_sig_n(i,j) = length(find(abs(fc_net_MedianWeight_n(i,j))<abs(fc_net_weight_per_n(:,i,j))))./size(fc_net_weight_per_n,1);
     end
 end

save([PLSca_Folder '/Weight_sig_Comp',num2str(component),'.mat'],'Behavior_Weight_actual', 'Behavior_weight_median','Behavior_weight_sig','Brain_Weight_actual', 'fc_roi_weight_actual','fc_roi_MedianWeight','Brain_roi_weight_sig','fc_net_weight_actual_sumn','fc_net_weight_actual_abs','fc_net_weight_actual_p','fc_net_weight_actual_n','fc_net_MedianWeight_sumn','fc_net_MedianWeight_abs','fc_net_MedianWeight_p','fc_net_MedianWeight_n','Brain_net_weight_sig_sumn','Brain_net_weight_sig_abs', 'Brain_net_weight_sig_p','Brain_net_weight_sig_n','-v7.3');



