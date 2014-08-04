clear all
cd C:\WC_Germany\

%gains in format [MP LF2 LF3 LF8].

%% add mec pyramidal cells to db
load pers_revised_dir
gains = [50 2e3 2e3 3e3;
    50 2e3 2e3 3e3;
    40 2e3 2e3 3e3;
    100 3e3 3e3 4e3;
    100 3e3 3e3 4e3;
    100 3e3 3e3 4e3;
    90 3e3 3e3 4e3;
    120 3e3 3e3 4e4;
    100 2e3 2e3 5e3;
    80 2e3 2e3 5e3;
    100 3e3 3e3 4e3;
    70 3e3 3e3 4e3;
    80 2e3 2e3 3e3;
    100 2e3 2e3 3e3;
    80 2e3 2e3 3e3;
    100 2e3 2e3 4e3;
    90 1e3 2e3 5e3;
    50 1e3 2e3 3e3;
    50 1e3 2e3 3e3;
    50 1e3 2e3 3e3;
    100 1e3 2e3 3e3];

for i = 1:17
    wc_db(i).directory = strcat('C:\wc_data',dir_array{i}(32:end));
end
for i = 18:21
    wc_db(i).directory = strcat('C:\wc_data',dir_array{i}(28:end));
end

for i = 1:21
    f_names{i} = f_names{i}(5:end);% get rid of "pyr" prefix
    wc_db(i).date = f_names{i};
    wc_db(i).celltype = 'pyr';
    wc_db(i).region = 'mec';
    wc_db(i).layer = 'l3';
    wc_db(i).drugs = 'none';
    wc_db(i).gains = gains(i,:);
end

%% add lec cells to db
gains = [100 2e3 2e3 2e3;
    50 1e3 2e3 2e3;
    80 2e3 2e3 2e3;
    100 1e3 2e3 1e3;
    50 1e3 2e3 3e3;
    50 1e3 2e3 3e3;
    50 1e3 2e3 3e3];

for i = 22:25
    dir_array{i} = strcat('C:\wc_data',dir_array{i}(29:end));
end
for i = 26:28
    dir_array{i} = strcat('C:\wc_data',dir_array{i}(28:end));
end

for i = 22:28
    dblen = length(wc_db);
    f_names{i} = f_names{i}(5:end);% get rid of "lec" prefix
    wc_db(dblen+1).date = f_names{i};
    wc_db(dblen+1).directory = dir_array{i};
    wc_db(dblen+1).celltype = 'pyr';
    wc_db(dblen+1).region = 'lec';
    wc_db(dblen+1).layer = 'l3';
    wc_db(i).drugs = 'none';
    wc_db(i).gains = gains(i-21,:);
end
wc_db(25).layer = 'un'; %set 2005-12-09_A to unkown layer

gains = [80 1e3 2e3 1e3;
    80 1e3 2e3 1e3;
    60 1e3 2e3 1e3;
    70 1e3 2e3 2e3;
    70 1e3 2e3 2e3;
    50 1e3 2e3 2e3;
    50 1e3 2e3 3e3];

%now load other lec cells
already_used = [1 4 10 11];
load('C:\WC_Germany\lateralEC\LEC_dir.mat')
dir_array(already_used) = [];
f_names(already_used) = [];

for i = 1:length(dir_array)
    dir_array{i} = strcat('C:\wc_data',dir_array{i}(29:end));
end
for i = 1:length(dir_array)
    dblen = length(wc_db);
    wc_db(dblen+1).directory = dir_array{i};
    wc_db(dblen+1).date = f_names{i};
    wc_db(dblen+1).region = 'lec';
    wc_db(i).drugs = 'none';
    wc_db(i).gains = gains(i,:);
end

cur_loc = find_in_struct_vec(wc_db,'date','2005-11-25_B');
wc_db(cur_loc).layer = 'l2';
wc_db(cur_loc).celltype = 'un';

cur_loc = find_in_struct_vec(wc_db,'date','2005-11-28');
wc_db(cur_loc).layer = 'l2';
wc_db(cur_loc).celltype = 'ste';

cur_loc = find_in_struct_vec(wc_db,'date','2005-12-09_B');
wc_db(cur_loc).layer = 'l2';
wc_db(cur_loc).celltype = 'ste';

cur_loc = find_in_struct_vec(wc_db,'date','2005-12-09_C');
wc_db(cur_loc).layer = 'un';
wc_db(cur_loc).celltype = 'un';

cur_loc = find_in_struct_vec(wc_db,'date','2005-12-09_D');
wc_db(cur_loc).layer = 'un';
wc_db(cur_loc).celltype = 'un';

cur_loc = find_in_struct_vec(wc_db,'date','2005-12-12_A');
wc_db(cur_loc).layer = 'l3';
wc_db(cur_loc).celltype = 'pyr';

cur_loc = find_in_struct_vec(wc_db,'date','2005-12-12_B');
wc_db(cur_loc).layer = 'un';
wc_db(cur_loc).celltype = 'un';

%% add adjacent cells to db
gains = [50 1e3 2e3 2e3;
    50 1e3 2e3 3e3;
    40 2e3 2e3 3e3;
    50 2e3 2e3 3e3;
    100 2e3 2e3 3e3;
    nan nan nan nan];

load pers_revised_dir
for i = 29:34
    dir_array{i} = strcat('C:\wc_data',dir_array{i}(40:end));
end
for i = 29:34
    f_names{i} = f_names{i}(4:end);
    dblen = length(wc_db);
    wc_db(dblen+1).directory = dir_array{i};
    wc_db(dblen+1).date = f_names{i};
    wc_db(dblen+1).celltype = 'pyr';
    wc_db(dblen+1).region = 'adj';
    wc_db(dblen+1).layer = 'l3';
    wc_db(dblen+1).drugs = 'none';
    wc_db(dblen+1).gains = gains(i-28,:);
end 



%% add stellates to db
load('C:\WC_Germany\JMM_analysis_ste\dir_tree_ste.mat')
gains = [100 1e3 2e3 3e3;
    100 1e3 2e3 3e3;
    100 3e3 3e3 4e3;
    100 3e3 3e3 4e3;
    150 3e3 3e3 4e3;
    100 3e3 3e3 4e3;
    150 3e3 3e3 4e3;
    120 3e3 3e3 4e3;
    150 3e3 3e3 5e3;
    120 3e3 3e3 5e3;
    100 3e3 3e3 4e3;
    100 3e3 3e3 4e3;
    80 3e3 3e3 4e3;
    100 2e3 2e3 5e3;
    110 2e3 2e3 3e3;
    100 1e3 1e3 5e3];

for i = 1:length(f_names)
    dir_array{i} = strcat('C:\wc_data',dir_array{i}(32:end));
end

for i = 1:length(f_names)
    f_names{i} = f_names{i}(7:end); %get rid of "stell" prefix
    dblen = length(wc_db);
    wc_db(dblen+1).directory = dir_array{i};
    wc_db(dblen+1).date = f_names{i};
    wc_db(dblen+1).celltype = 'ste';
    wc_db(dblen+1).region = 'mec';
    wc_db(dblen+1).layer = 'l2';
    wc_db(i).drugs = 'none';
    wc_db(i).gains = gains(i,:);
end

%% add layer 5 cells to db
load('C:\WC_Germany\Layer5\layer_5_dir.mat')
gains = [100 2e3 2e3 5e3;
    100 2e3 2e3 5e3;
    100 1e3 2e3 5e3;
    100 1e3 2e3 5e3;
    100 1e3 2e3 5e3;
    100 1 1 5e3;
    100 1 1e3 5e3;
    90 1 1e3 5e3;
    100 1 1e3 5e3;
    100 1 1e3 5e3;
    80 1e3 1e3 5e3;
    100 1e3 1e3 5e3;
    100 1e3 1e3 5e3;
    100 1e3 1e3 5e3;
    100 1e3 1e3 5e3;
    90 1e3 1e3 5e3;
    90 1e3 1e3 5e3];

for i = 1:length(f_names)
    dir_array{i} = strcat('C:\wc_data',dir_array{i}(41:end));
end

for i = 1:length(f_names)
    dblen = length(wc_db);
    wc_db(dblen+1).directory = dir_array{i};
    wc_db(dblen+1).date = f_names{i};
    wc_db(dblen+1).celltype = 'pyr';
    wc_db(dblen+1).region = 'mec';
    wc_db(dblen+1).layer = 'l5';
    wc_db(i).drugs = 'none';
    wc_db(i).gains = gains(i,:);
end

%% add atropine data to db

load('C:\WC_Germany\Atropine\atropine_dir.mat')
gains = [nan nan nan nan;
     nan nan nan nan;
    nan nan nan nan;
    nan nan nan nan;
    nan nan nan nan;
    nan nan nan nan;
    nan nan nan nan;
    nan nan nan nan;
    nan nan nan nan];

atr_type = [7 7 8 8 8 7 7 7 8];

for i = 1:length(f_names)
    dir_array{i} = strcat('C:\wc_data',dir_array{i}(21:end));
end
for i = 1:length(f_names)
    dblen = length(wc_db);
    wc_db(dblen+1).directory = dir_array{i};
    wc_db(dblen+1).date = f_names{i};
    wc_db(dblen+1).celltype = 'un';
    wc_db(dblen+1).region = 'mec';
    wc_db(dblen+1).layer = 'un';
    wc_db(i).drugs = 'atropine';
    wc_db(i).gains = gains(i,:);
end

cur_loc = find_in_struct_vec(wc_db,'date','2007-11-22_A');
wc_db(cur_loc).layer = 'l3';
wc_db(cur_loc).celltype = 'pyr';

cur_loc = find_in_struct_vec(wc_db,'date','2007-11-22_B');
wc_db(cur_loc).layer = 'l3';
wc_db(cur_loc).celltype = 'pyr';

cur_loc = find_in_struct_vec(wc_db,'date','2007-11-26_B');
wc_db(cur_loc).layer = 'l2';
wc_db(cur_loc).celltype = 'ste';

cur_loc = find_in_struct_vec(wc_db,'date','2007-11-27_A');
wc_db(cur_loc).layer = 'l2';
wc_db(cur_loc).celltype = 'ste';

cur_loc = find_in_struct_vec(wc_db,'date','2007-11-27_B');
wc_db(cur_loc).layer = 'l2';
wc_db(cur_loc).celltype = 'ste';

cur_loc = find_in_struct_vec(wc_db,'date','2007-11-27_C');
wc_db(cur_loc).layer = 'l3';
wc_db(cur_loc).celltype = 'pyr';

cur_loc = find_in_struct_vec(wc_db,'date','2007-11-27_D');
wc_db(cur_loc).layer = 'l3';
wc_db(cur_loc).celltype = 'pyr';

cur_loc = find_in_struct_vec(wc_db,'date','2007-11-28_A');
wc_db(cur_loc).layer = 'l3';
wc_db(cur_loc).celltype = 'pyr';

cur_loc = find_in_struct_vec(wc_db,'date','2007-11-28_B');
wc_db(cur_loc).layer = 'l2';
wc_db(cur_loc).celltype = 'ste';


%% add parietal cortical cells to db
load('C:\WC_Germany\overall_calcs\HC_Cort\hc_cort_dir_temp')
num_cort_cells = 11;
gains = [70 1e3 2e3 1e3;
    60 2e3 5e3 2e3;
    100 1e3 2e3 1e3;
    80 1e3 2e3 1e3;
    60 1e3 2e3 1e3;
    50 1e3 2e3 2e3;
    70 1e3 2e3 2e3;
    50 1e3 2e3 2e3;
    50 1e3 2e3 2e3;
    60 1e3 2e3 2e3;
    60 1e3 2e3 2e3];

for i = 1:1
    dir_array{i} = strcat('C:\wc_data',dir_array{i}(19:end));
end
for i = 2:2
    dir_array{i} = strcat('C:\wc_data',dir_array{i}(43:end));
end
for i = 3:num_cort_cells
    dir_array{i} = strcat('C:\wc_data',dir_array{i}(19:end));
end
for i = 1:num_cort_cells
    dblen = length(wc_db);
    wc_db(dblen+1).directory = dir_array{i};
    wc_db(dblen+1).date = f_names{i};
    wc_db(dblen+1).celltype = 'pyr';
    wc_db(dblen+1).region = 'par';
    wc_db(dblen+1).layer = 'l23';
    wc_db(i).drugs = 'none';
    wc_db(i).gains = gains(i,:);
end
dir_array(1:num_cort_cells) = [];
f_names(1:num_cort_cells) = [];

%% add ca1 pyr cells to db
num_ca1_cells = 14;
gains = [100 2e3 5e3 2e3;
    nan nan 8e3 2.5e3;
    100 2e3 6e3 2e3;
    100 2e3 2e3 2e3;
    200 1e3 2e3 2e3;
    200 1e3 1e3 1e3;
    100 1e3 1e3 1e3;
    80 1e3 2e3 2e3;
    80 1e3 2e3 2e3;
    80 1e3 2e3 2e3;
    80 2e3 2e3 2e3;
    50 1e3 2e3 1e3;
    80 1e3 2e3 3e3;
    50 1e3 2e3 3e3];
for i = 1:4
    dir_array{i} = strcat('C:\wc_data',dir_array{i}(27:end));
end
for i = 5:num_ca1_cells
    dir_array{i} = strcat('C:\wc_data',dir_array{i}(19:end));
end

for i = 1:num_ca1_cells
    dblen = length(wc_db);
    wc_db(dblen+1).directory = dir_array{i};
    wc_db(dblen+1).date = f_names{i};
    wc_db(dblen+1).celltype = 'pyr';
    wc_db(dblen+1).region = 'ca1';
    wc_db(dblen+1).layer = 'na';
    wc_db(i).drugs = 'none'; 
    wc_db(i).gains = gains(i,:);
end
dir_array(1:num_ca1_cells) = [];
f_names(1:num_ca1_cells) = [];

%% add ca1 int cells to db
num_ca1_int_cells = 8;
gains = [70 1e3 2e3 1e3;
    60 3e3 5e3 3e3;
    80 1e3 1e3 1e3;
    80 1e3 2e3 2e3;
    100 2e3 2e3 3e3;
    60 1e3 2e3 1e3;
    nan nan nan nan;
    100 1e3 2e3 3e3];

for i = 1:1
    dir_array{i} = strcat('C:\wc_data',dir_array{i}(19:end));
end
for i = 2:2
    dir_array{i} = strcat('C:\wc_data',dir_array{i}(27:end));
end
for i = 3:num_ca1_int_cells
    dir_array{i} = strcat('C:\wc_data',dir_array{i}(19:end));
end

for i = 1:num_ca1_int_cells
    dblen = length(wc_db);
    wc_db(dblen+1).directory = dir_array{i};
    wc_db(dblen+1).date = f_names{i};
    wc_db(dblen+1).celltype = 'int';
    wc_db(dblen+1).region = 'ca1';
    wc_db(dblen+1).layer = 'na';
    wc_db(i).drugs = 'none'; 
    wc_db(i).gains = gains(i,:);
end
dir_array(1:num_ca1_int_cells) = [];
f_names(1:num_ca1_int_cells) = [];

%% add ca3 pyr cells to db
num_ca3_cells = 10;
gains = [100 1e3 2e3 3e3;
    100 1e3 2e3 2e3;
    100 1e3 2e3 4e3;
    100 2e3 2e3 3e3;
    100 1e3 1e3 1e3;
    150 2e3 2e3 3e3;
    100 2e3 2e3 3e3;
    90 2e3 2e3 3e3;
    100 3e3 3e3 4e3;
    150 2e3 3e3 4e3];

for i = 1:4
    dir_array{i} = strcat('C:\wc_data',dir_array{i}(30:end));
end
for i = 5:num_ca3_cells
    dir_array{i} = strcat('C:\wc_data',dir_array{i}(19:end));
end

for i = 1:num_ca3_cells
    dblen = length(wc_db);
    wc_db(dblen+1).directory = dir_array{i};
    wc_db(dblen+1).date = f_names{i};
    wc_db(dblen+1).celltype = 'pyr';
    wc_db(dblen+1).region = 'ca3';
    wc_db(dblen+1).layer = 'na';
    wc_db(i).drugs = 'none'; 
    wc_db(i).gains = gains(i,:);
end
dir_array(1:num_ca3_cells) = [];
f_names(1:num_ca3_cells) = [];

%% add dg cells to db
num_dg_cells = 10;
gains = [100 1e3 2e3 2e3;
    50 1e3 2e3 3e3;
    100 2e3 2e3 3e3;
    100 2e3 2e3 3e3;
    nan nan nan nan;
    100 1e3 2e3 3e3;
    100 2e3 2e3 3e3;
    150 2e3 2e3 3e3;
    200 3e3 3e3 4e3;
    300 2e3 3e3 4e4];
  
for i = 1:1
    dir_array{i} = strcat('C:\wc_data',dir_array{i}(30:end));
end
for i = 2:2
    dir_array{i} = strcat('C:\wc_data',dir_array{i}(29:end));
end
for i = 3:num_dg_cells
    dir_array{i} = strcat('C:\wc_data',dir_array{i}(19:end));
end

for i = 1:num_dg_cells
    dblen = length(wc_db);
    wc_db(dblen+1).directory = dir_array{i};
    wc_db(dblen+1).date = f_names{i};
    wc_db(dblen+1).celltype = 'gran';
    wc_db(dblen+1).region = 'dg';
    wc_db(dblen+1).layer = 'na';
    wc_db(i).drugs = 'none'; 
end


%% now save
cd C:\WC_Germany
save wc_database wc_db