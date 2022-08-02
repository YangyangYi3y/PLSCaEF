clear
Folderdir = 'D:\projects\ABCD\func_baseline_10min\ABCDfunc10minbaseline';
Foldername = dir(Folderdir);
for i = 3:length(Foldername)
    fc_dir = [Folderdir,'\',Foldername(i).name];
    subjectkey1 = split(fc_dir,'-');
    subjectkey2 = split(subjectkey1{2},'_');
    subjectkeyi = subjectkey2 {1};
    fc1 = cifti_read(fc_dir);
    fc2 = fc1.cdata;
    fc3=fc2-diag(diag(fc2));
    fci=fc2([triu(fc2)-diag(diag(fc2))]~=0)';
    subjectkey{i-2} = subjectkeyi ;
    fc(i-2,:)=(fci);
end

subjectkey = subjectkey';
fct = table(subjectkey,fc,'VariableNames',{'subjectkey','fc'});
save('D:\projects\ABCD\func_baseline_10min\participants_name.txt','subjectkey');
save('D:\projects\ABCD\func_baseline_10min\fc.mat','fc');