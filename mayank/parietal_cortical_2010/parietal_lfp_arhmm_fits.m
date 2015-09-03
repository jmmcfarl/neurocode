%inputs:
% raw_data
% raw_Fs
clear all
close all

drive_letter = 'F';

addpath(strcat(drive_letter,':\Code\smoothing\software'))
addpath(strcat(drive_letter,':\Code\FullBNT-1.0.4\KPMstats\'))
addpath(strcat(drive_letter,':\Code\FullBNT-1.0.4\netlab3.3'))
addpath(strcat(drive_letter,':\WC_Germany\hsmm_state_detection'))
addpath(strcat(drive_letter,':\WC_Germany\parietal_cortical_2010\'))
addpath(strcat(drive_letter,':\Code\maphmmbox\'))
addpath(strcat(drive_letter,':\Code\arfit'))
cd(strcat(drive_letter,':\WC_Germany\parietal_cortical_2010\'))
load parietal_cortical_2010
load desynch_times_individual
raw_Fs = 2016;

% get rid of interneurons
interneurons = find_struct_field_vals(sess_data,'cell_type','interneuron');
sess_data(interneurons) = [];
desynch_times_mp(interneurons) = [];
desynch_times_lf8(interneurons) = [];
desynch_times_lf4(interneurons) = [];


frontal = find_struct_field_vals(sess_data,'region','frontal');
prefrontal = find_struct_field_vals(sess_data,'region','prefrontal');
parietal = find_struct_field_vals(sess_data,'region','parietal');
thom_el = find_struct_field_vals(sess_data,'thom_elec',1);
thom_par = thom_el(find(ismember(thom_el,parietal)));
thom_pfc = setdiff(thom_el,thom_par);

sess_data = sess_data(thom_el);
desynch_times_lf4 = desynch_times_lf4(thom_el);

n = length(sess_data);
for d = 1:n
    fprintf('session %d\n',d)
    cdir = sess_data(d).directory;
    cdir(1) = drive_letter;
    cd(cdir)
    load used_data wcv_minus_spike lf4

   [armm_state_seq4,arhmm4,Fs_bb,Fs,gamma] = ...
            parietal_get_arhmm_uds_state_seq(lf4,raw_Fs,desynch_times_lf4{d},strcat(sess_data(d).name,'_lf4'));

end

cd F:\WC_Germany\parietal_cortical_2010
