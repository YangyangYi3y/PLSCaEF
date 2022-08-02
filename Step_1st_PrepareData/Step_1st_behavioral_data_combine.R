rm(list=ls())
datapath <- '/GPFS/cuizaixu_lab_permanent/yiyangyang/projects/ABCD/check_data/'

#1 read tbss 
tbbs <- read.table(paste(datapath,'abcd_tbss01','.txt',sep=''),
                   stringsAsFactors = FALSE,na.strings="",header = TRUE)
tbbs_baseline <- tbbs[tbbs$eventname =='baseline_year_1_arm_1',]
tbbs_variable_names<-c('subjectkey','nihtbx_flanker_uncorrected',
                  'nihtbx_list_uncorrected',
                  'nihtbx_cardsort_uncorrected',
                  'nihtbx_pattern_uncorrected',
                  'nihtbx_picture_uncorrected')

tbbs_data <- tbbs_baseline[tbbs_variable_names]


#2 processing lmtp201
lmtp201 <- read.table(paste(datapath,'lmtp201','.txt',sep=''),
                      stringsAsFactors = FALSE,na.strings="",header = TRUE)
lmtp201_baseline <- lmtp201[lmtp201$eventname =='baseline_year_1_arm_1',]
lmtp201_names<-c('subjectkey','lmt_scr_perc_correct','lmt_scr_rt_correct')
lmtp201_datapre1 <- lmtp201_baseline[lmtp201_names]
lmtp201_datapre2 <- na.omit(lmtp201_datapre1)
lmtp201_datapre3<- transform(lmtp201_datapre2,
                             lmt_scr_efficiency = 
                             as.numeric(lmtp201_datapre2$lmt_scr_perc_correct)/as.numeric(lmtp201_datapre2$lmt_scr_rt_correct))

lmtp201_variable_names<-c('subjectkey','lmt_scr_efficiency')
lmtp201_data <- lmtp201_datapre3[lmtp201_variable_names]

#3 read cct01
cct <- read.table(paste(datapath,'cct01','.txt',sep=''),
                   stringsAsFactors = FALSE,na.strings="",header = TRUE)
cct_baseline <- cct[cct$eventname =='baseline_year_1_arm_1',]
cct_variable_names<-c('subjectkey','cash_choice_task')
cct_data <- cct_baseline[cct_variable_names]


#4 processing abcd_ps01
abcd_ps01 <- read.table(paste(datapath,'abcd_ps01','.txt',sep=''),
                  stringsAsFactors = FALSE,na.strings="",header = TRUE)
abcd_ps01_baseline <- abcd_ps01[abcd_ps01$eventname =='baseline_year_1_arm_1',]
ps_variable_names<-c('subjectkey','pea_wiscv_trs',
                     'pea_ravlt_sd_trial_i_tc',
                     'pea_ravlt_sd_trial_ii_tc',
                     'pea_ravlt_sd_trial_iii_tc',
                     'pea_ravlt_sd_trial_iv_tc',
                     'pea_ravlt_sd_trial_v_tc',
                     'pea_ravlt_sd_trial_vi_tc',
                     'pea_ravlt_ld_trial_vii_tc',
                     'pea_ravlt_sd_listb_tc')
ps <- abcd_ps01_baseline[ps_variable_names]
ps <- na.omit(ps)
ps_num <- apply(ps[,2:10],2,as.numeric)
ps_num$pea_ravlt_sd_listb_tc
pea_ravlt_sd_listb_tc
ps <- transform(ps,
                pi = as.numeric(ps$pea_ravlt_sd_listb_tc)/as.numeric(ps$pea_ravlt_sd_trial_i_tc),
                ri = as.numeric(ps$pea_ravlt_sd_trial_vi_tc)/as.numeric(ps$pea_ravlt_sd_trial_v_tc),
                lr = (as.numeric(ps$pea_ravlt_sd_trial_i_tc)+
                        as.numeric(ps$pea_ravlt_sd_trial_ii_tc)+
                        as.numeric(ps$pea_ravlt_sd_trial_iii_tc)+
                        as.numeric(ps$pea_ravlt_sd_trial_iv_tc)+
                        as.numeric(ps$pea_ravlt_sd_trial_v_tc))
                      -5*as.numeric(ps$pea_ravlt_sd_trial_i_tc),
                fs = as.numeric(ps$pea_ravlt_ld_trial_vii_tc)/as.numeric(ps$pea_ravlt_sd_trial_vi_tc),
                mr = (as.numeric(ps$pea_ravlt_sd_trial_i_tc)+
                      as.numeric(ps$pea_ravlt_sd_trial_ii_tc)+
                      as.numeric(ps$pea_ravlt_sd_trial_iii_tc)+
                      as.numeric(ps$pea_ravlt_sd_trial_iv_tc)+
                      as.numeric(ps$pea_ravlt_sd_trial_v_tc))-35)
abcd_ps01_variable_names<-c('subjectkey','pi','ri','lr','fs','mr')

abcd_ps01_data <- ps[abcd_ps01_variable_names]



#5 processing mribrec02
mribrec02 <- read.table(paste(datapath,'mribrec02','.txt',sep=''),
                      stringsAsFactors = FALSE,na.strings="",header = TRUE)
mribrec02_baseline <- mribrec02[mribrec02$eventname =='baseline_year_1_arm_1',]
mribrec02_variable_names<-c('subjectkey'
                            ,'tfmri_rec_all_beh_posf_dpr',
                            'tfmri_rec_all_beh_neutf_dp',
                            'tfmri_rec_all_beh_negf_dp',
                            'tfmri_rec_all_beh_place_dp')
mribrec02_data <- mribrec02_baseline[mribrec02_variable_names]


# preprocession mri task data, add inclusion criteria
inclusion_criteria <- read.table(paste(datapath,'abcd_imgincl01.txt',sep=''),
                                 stringsAsFactors = FALSE,na.strings="",header = TRUE)
abcd_mrinback02_incri <- inclusion_criteria[inclusion_criteria$imgincl_nback_include==1,]
abcd_mrinback02_ob <- abcd_mrinback02_incri$subjectkey
abcd_sst02_incri <- inclusion_criteria[inclusion_criteria$imgincl_sst_include==1,]
abcd_sst02_ob <- abcd_sst02_incri$subjectkey
abcd_mid02_incri <- inclusion_criteria[inclusion_criteria$imgincl_mid_include ==1,]
abcd_mid02_ob <- abcd_mid02_incri$subjectkey

#6 processing abcd_sst02
abcd_sst02 <- read.csv(paste(datapath,'abcd_sst02','.csv',sep=''),
                      stringsAsFactors = FALSE,na.strings="",header = TRUE)
abcd_sst02_baseline <- abcd_sst02 [abcd_sst02$eventname =='baseline_year_1_arm_1',]
abcd_sst02_baseline_inclu <- abcd_sst02_baseline[abcd_sst02_baseline$subjectkey %in% abcd_sst02_ob,]
abcd_sst02_variable_names<-c('subjectkey','tfmri_sst_all_beh_total_mssrt')
abcd_sst02_data <- abcd_sst02_baseline_inclu[abcd_sst02_variable_names]


#7 processing abcd_mrinback02
abcd_mrinback02 <- read.table(paste(datapath,'abcd_mrinback02','.txt',sep=''),
                       stringsAsFactors = FALSE,na.strings="",header = TRUE)
abcd_mrinback02_baseline <- abcd_mrinback02[abcd_mrinback02$eventname =='baseline_year_1_arm_1',]
abcd_mrinback02_baseline_inclu <- abcd_mrinback02_baseline[abcd_mrinback02_baseline$subjectkey 
                                                           %in% abcd_mrinback02_ob,]
abcd_mrinback02_variable_names<-c('subjectkey','tfmri_nb_all_beh_c2b_rate','tfmri_nb_all_beh_c2b_mrt')
abcd_mrinback02_data <- abcd_mrinback02_baseline_inclu[abcd_mrinback02_variable_names]

#8 processsion abcd_mid02
abcd_mid02 <- read.table(paste(datapath,'abcd_mid02','.txt',sep=''),
                              stringsAsFactors = FALSE,na.strings="",header = TRUE)
abcd_mid02_baseline <- abcd_mid02[abcd_mid02$eventname =='baseline_year_1_arm_1',]
abcd_mid02_baseline_inclu <- abcd_mid02_baseline[abcd_mid02_baseline$subjectkey 
                                                           %in% abcd_mid02_ob,]
abcd_mid02_variable_names<-c('subjectkey','tfmri_mid_all_beh_srwpfb_mrt',
                             'tfmri_mid_all_beh_lrwpfb_mrt',
                             'tfmri_mid_all_beh_slpfb_mrt',
                             'tfmri_mid_all_beh_llpfb_mrt')
abcd_mid02_data <- abcd_mid02_baseline_inclu[abcd_mid02_variable_names]   

#merge data
ef <- merge(tbbs_data,lmtp201_data, by = 'subjectkey')
ef <- merge(ef, cct_data, by = 'subjectkey')
ef <-merge(ef, abcd_ps01_data, by = 'subjectkey')
ef <-merge(ef, mribrec02_data, by = 'subjectkey')
ef <-merge(ef, abcd_mrinback02_data,by = 'subjectkey')
ef <-merge(ef, abcd_sst02_data, by = 'subjectkey')
ef <-merge(ef, abcd_mid02_data,by = 'subjectkey')

ef_martrix = as.matrix(ef)
ef_num <- apply(ef_martrix[,2:24],2,as.numeric)
index = rowSums(is.finite(ef_num))
ef_inf = ef[index == 23,]
for (i in nrow(ef_inf)){
  ef_subjectkey<- gsub('[_]','',ef_inf[,1])
  ef_inf$subjectkey <-ef_subjectkey
}
ef_num_inf <- ef_num[index ==23,]

subjectkey <- read.table('D:/projects/ABCD/func_baseline_10min/subjectkey.txt',
                       stringsAsFactors = FALSE,na.strings="",header = TRUE)

mri_subjectkey<- data.frame(subjectkey)
ef_z_inf <- scale(ef_num_inf)
#repalce outlier with 3 or -3
ef_z_inf[ef_z_inf > 3] = 3 
ef_z_inf[ef_z_inf < -3] = -3 
ef_zdata = as.data.frame(ef_z_inf)
ef_zdata$subjectkey <- ef_inf$subjectkey
efz_withfmri <-merge(ef_zdata,mri_subjectkey,by = 'subjectkey')   
write.csv(efz_withfmri, file = "D:/projects/ABCD/EF/data/ef_task23.csv")


