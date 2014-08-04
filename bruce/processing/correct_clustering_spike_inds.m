clear all

global data_dir base_save_dir init_save_dir Expt_name Vloaded n_probes loadedData
Expt_name = 'M232';

% Expt_name = 'G086';
% data_loc = '/media/NTlab_data1/Data/bruce/';
data_loc = '/home/james/Data/bruce/';

base_save_dir = ['~/Analysis/bruce/' Expt_name '/clustering'];
init_save_dir = ['~/Analysis/bruce/' Expt_name '/clustering/init'];

if ~exist(base_save_dir,'dir');
    mkdir(base_save_dir);
end
if ~exist(init_save_dir,'dir');
    mkdir(init_save_dir);
end

%location of FullV files
data_dir = [data_loc Expt_name];

%location of Expts.mat files
data_dir2 = ['~/Data/bruce/' Expt_name];
% data_dir2 = ['/media/NTlab_data1/Data/bruce/' Expt_name];

Vloaded = nan;
n_probes = 24;

cd(data_dir2);
if Expt_name(1) == 'G';
    if strcmp(Expt_name,'G029')
        load('G029Expts.mat');
    else
        load(sprintf('jbe%sExpts.mat',Expt_name));
    end
    n_probes = 96;
elseif Expt_name(1) == 'M'
    load(sprintf('lem%sExpts.mat',Expt_name));
    n_probes = 24;
end
n_blocks = length(Expts);
target_probes = 1:n_probes;
%%
for ii = 1:n_blocks
    %load existing clusters for this block
    cur_dat_name = [base_save_dir sprintf('/Block%d_Clusters.mat',ii)];
    clear loadedData
    if exist(cur_dat_name,'file')
        %for LP load all Voltage signals for this block
        if Expt_name(1) == 'M'
            sfile_name = [data_dir sprintf('/Expt%dFullV.mat',ii)];
            if Vloaded ~= ii
                fprintf('Loading data file %s\n',sfile_name);
                [loadedData.V,loadedData.Vtime,loadedData.Fs] = Load_FullV(sfile_name, false, [100 nan],1:n_probes);
                Vloaded = ii;
            end
        elseif Expt_name(1) == 'G'
            sfile_name = [data_dir sprintf('/Expt%d.p%dFullV.mat',ii,1)];
                [loadedData.V,loadedData.Vtime,loadedData.Fs] = Load_FullV(sfile_name, false, [100 nan]);
        end

        
        fprintf('Loading clusters for block %d\n',ii);
        load(cur_dat_name);
       
        for probe_num = target_probes
            cur_spk_times = Clusters{probe_num}.times;
            cur_spk_inds = Clusters{probe_num}.spk_inds;
%             new_spk_inds2 = round(interp1(loadedData.Vtime,1:length(loadedData.Vtime),cur_spk_times));
            new_spk_inds = find(ismember(loadedData.Vtime,cur_spk_times));
            Clusters{probe_num}.spk_inds = new_spk_inds;
            if length(new_spk_inds) ~= length(cur_spk_times)
                error('mismatch')
            end
            if length(cur_spk_times) ~= length(new_spk_inds) | ...
                    length(cur_spk_times) ~= length(Clusters{probe_num}.spike_clusts)
                error('mismatch')
            end
            dspk = setdiff(cur_spk_inds,new_spk_inds);
%             dspk2 = setdiff(new_spk_inds,new_spk_inds2);
            fprintf('Probe %d, %d of %d diff spikes\n',probe_num,length(dspk),length(cur_spk_times));
%             fprintf('Probe %d, %d of %d diff2 spikes\n',probe_num,length(dspk2),length(cur_spk_times));
        end
%         save(cur_dat_name,'Clusters');
        clear Clusters
    end
    
end
