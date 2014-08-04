
%% create master data directory
% cell_type variable organized as follows:
% 1 = MEC L3 pyr
% 2 = MEC L2 ste
% 3 = L5
% 4 = LEC L3 pyr
% 5 = LEC L2 ste
% 6 = LEC unknown
% 7 = MEC L3 pyr atr
% 8 = MEC L2 ste atr

clear all
close all

%% load L3 pyramidal directory
load('C:\WC_Germany\JMM_analysis_pyr\dir_tree_update.mat')
over_dir = dir_array;
over_names = f_names;
cell_type = ones(1,length(dir_array));
clear dir_array f_names

%% load L2 stellate directory
load('C:\WC_Germany\JMM_analysis_ste\dir_tree_ste.mat')
over_dir = [over_dir dir_array];
over_names = [over_names f_names];
cell_type = [cell_type 2*ones(1,length(dir_array))];
clear dir_array f_names

%% load L5 directory
load('C:\WC_Germany\Layer5\layer_5_dir.mat')
over_dir = [over_dir dir_array];
over_names = [over_names f_names];
cell_type = [cell_type 3*ones(1,length(dir_array))];
clear dir_array f_names

%% load LEC directory
lec_type = [4 6 5 5 5 6 6 4 6 4 4];
load('C:\WC_Germany\lateralEC\LEC_dir.mat')
over_dir = [over_dir dir_array];
over_names = [over_names f_names];
cell_type = [cell_type lec_type];
clear dir_array f_names

%% load atropine directory
load('C:\WC_Germany\Atropine\atropine_dir.mat')
atr_type = [7 7 8 8 8 7 7 7 8];
over_dir = [over_dir dir_array];
over_names = [over_names f_names];
cell_type = [cell_type atr_type];
clear dir_array f_names

%%
cd C:\WC_Germany\overall_calcs
clear d ans
save overall_dir