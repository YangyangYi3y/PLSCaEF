clear
%get fd data
%calculate mean framewise displacement
fold_path = 'D:\projects\ABCD\func_baseline_10min\abcd_FD';
file_name = [fold_path,'/desc-filteredincludingFD_motion'];
%get subjects' name
file_namedir = dir(file_name);
for i = 3:length(file_namedir)
    fd_name{i-2} = file_namedir(i).name;
    sub_namepre1 {i-2}= file_namedir(i).name;
    sub_namepre2{i-2} = split(sub_namepre1 {i-2},'-');
    sub_namepre3{i-2}  =  split(sub_namepre2{i-2}(2),'_');
    sub_name(i-2) = sub_namepre3{i-2}(1);   
end
sub_name= sub_name';
sub_name_unique = unique(sub_name);
%read FD included file
for i = 1:length(fd_name)
    fd_namei = fd_name(i);
    [a] = textread([file_name,'\',fd_namei{1}],'%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%s');
    fd = str2double(a(3:end));
    fd_sum = sum(fd);
    fd_number = length(fd);
    sum_fd_total(i,:) = [sub_name(i),fd_sum,fd_number]; 
end

%calulate mean FD for each participante
for i = 1:length(sub_name_unique)
   a = sub_name_unique(i);
   index_name = strmatch(a,sub_name);
   sum_fdipre = sum_fd_total(index_name',:);
   sum_fdi = sum(cell2mat(sum_fdipre(:,2)));
   num_fdipre = sum(cell2mat(sum_fdipre(:,3)));
   mean_fdi = sum_fdi/num_fdipre;
   mean_fd(i) = mean_fdi; 
end
mean_fd = mean_fd';  
FD_subname= sub_name_unique;

% merge age sex fd batch by subjectkey
a = importdata('D:\projects\ABCD\participants_v1.0.0ORIGINAL\participant_info.txt');
participnat_id = a.textdata(2:end,1);
participnat_id = split(participnat_id,'-');
participnat_id = participnat_id(:,2);
site =  a.textdata(2:end,4);
scanner =  a.textdata(2:end,5);
age = a.data(:,4);
sex = a.data(:,2);
%888 is missing data, delete missing data
id = [find(strcmp(site, '888' ));find(strcmp(scanner, '888' ))];
participnat_id(id,:) = [];
age(id,:) = [];
sex(id,:) = [];
scanner(id,:) = [];
site(id,:) = [];

for i = 1:length(site)
    site_scanner{i} = [site{i},scanner{i}];
end
site_scanner_unique = unique(site_scanner);

batch = {};
for i = 1:length(site_scanner)
    batch(i,1) = cellstr(participnat_id{i});
    batch(i,2) = num2cell(find(strcmp(site_scanner_unique,site_scanner{i})));
    batch(i,3) = cellstr(site{i});
    batch(i,4)= cellstr(scanner{i});
end


% merge age sex fd by subjectkey
ef = importdata(['D:\projects\ABCD\EF\data\ef_task23.csv']);
EF_Data = ef.data;
EF_ItemName = ef.textdata(1,3:end);
EF_SubjectKey = ef.textdata(2:end, 2);

inter_subkey = intersect(EF_SubjectKey,participnat_id);
inter_subkey = intersect(FD_subname,inter_subkey)
cov = {};
for i = 1:length(inter_subkey)
  ind = find(strcmp(participnat_id,(inter_subkey{i})));
  %some subjects are replicated in the participnat_id, only extract the first one
  ind = ind(1);
  cov(i, 1) = inter_subkey(i, :);
  cov(i, 2) = num2cell(age(ind));
  cov(i, 3) = num2cell(sex(ind));
  cov(i, 4) = num2cell(batch{ind, 2});
end

inter_subkeyfd = intersect(FD_subname,inter_subkey);
for i = 1:length(inter_subkeyfd)
  ind = find(strcmp(FD_subname,(inter_subkeyfd{i})));
  %some subjects are replicated in the participnat_id, only extract the first one
  ind = ind(1);
  fd_subname(i) = FD_subname(ind);
  fd(i) = mean_fd(ind);
  ind = find(strcmp(FD_subname,(inter_subkeyfd{i})));
end
fd = fd';

age = cov(:, 2);
sex =cov(:, 3);
batch = cov(:, 4);

covariable = table(inter_subkeyfd,age,sex,fd,batch);
writetable(covariable,'D:\projects\ABCD\func_baseline_10min\covariables.txt')

