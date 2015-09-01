clear all
close all

addpath('F:\WC_Germany\parietal_cortical_2010\')
addpath('F:\WC_Germany\hsmm_state_detection\')

cd F:\WC_Germany\parietal_cortical_2010\
load parietal_cortical_2010
load F:\WC_Germany\parietal_cortical_2010\desynch_times_individual


frontal = find_struct_field_vals(sess_data,'region','frontal');
prefrontal = find_struct_field_vals(sess_data,'region','prefrontal');
parietal = find_struct_field_vals(sess_data,'region','parietal');

% sess_data = sess_data(parietal);
% desynch_times = desynch_times(parietal);

%get rid of interneurons
interneurons = find_struct_field_vals(sess_data,'cell_type','interneuron');
sess_data(interneurons) = [];
desynch_times_lf8(interneurons) = [];

raw_Fs = 2016;
dsf = 8;
Fsd = raw_Fs/dsf;

% thomas_el = find_struct_field_vals(sess_data,'thom_elec',1);
% sess_data = sess_data(thomas_el);
% desynch_times = desynch_times(thomas_el);

n = length(sess_data);

for d = 1:n
    
       cdir = sess_data(d).directory;
    cdir(1) = 'F';
    cd(cdir)
    pwd

    load hsmm_state_seq8_seg_lf_indivdesynch
    up_meanfun = [];
    down_meanfun = [];
    for i = 1:hmm8.Nsegs
        up_meanfun = [up_meanfun; hmm8.state(2).meanfun{i}];
        down_meanfun = [down_meanfun; hmm8.state(1).meanfun{i}];
    end
    temp = corrcoef(up_meanfun,down_meanfun);
    mean_corr(d) = temp(2,1);
end