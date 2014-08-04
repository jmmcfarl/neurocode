% ori = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170];
% sfreqCycDeg = [0.1, 0.15, 0.20, 0.25, 0.35, 0.5, 0.75, 1, 1.25, 1.5];
% phase0 = [0, 90, 180, 270];


cd ~/Data/blanche/rec_73/
clear all; load stims73.mat; load spks73.mat;

npix=16; tstids = stids20; tspks =spks20;  %stas for stim40 and stim20 are VERY similar

lutable = zeros(721,npix*npix);
for irow = 1:720;
    tor = ori(ctab(irow,2));
    tsf = 10*sfreqCycDeg(ctab(irow,3));
    tpha = phase0(ctab(irow,4));
    lutable(irow,:) = flatten(hartley3p(tor,tsf,tpha,npix));
end;

istim = zeros(length(tstids),npix^2);
for itime = 1:length(tstids);
    istim(itime,:) = lutable(tstids(itime)+1,:);
end;

%%
n_stims = 720;
stim_id_mat = cell(n_stims);
for j = 1:n_stims
    stim_id_mat{j} = find(tstids==j);
    n_stim_vec(j) = length(stim_id_mat{j});
end

%% CREATE TIME_EMBEDDED STIM MATRIX. NOT REALLY NECESSARY...
n_lags = 10;
tlags = (9:-1:0)*dt;
S = makeStimRows(tstids, n_lags, 0);
S(S==0) = 721;
%% CREATE FULL CELL RESPONSE MATRIX
n_cells = length(tspks);
cell_response_mat = zeros(n_cells,n_stims+1,n_lags);
for cc = 1:n_cells
    cc
    cur_spks = ceil(tspks{cc}/1000/dt);
    n_spks(cc) = length(cur_spks);
    for ss = 1:n_spks(cc)
        for nn = 1:n_lags
            cell_response_mat(cc,S(cur_spks(ss),nn),nn) = cell_response_mat(cc,S(cur_spks(ss),nn),nn) + 1;
        end
    end
end
cell_response_array = nan(n_cells,length(sfreqCycDeg),length(ori),length(phase0),n_lags);
for i = 1:720
    cell_response_array(:,ctab(i,3),ctab(i,2),ctab(i,4),:) = cell_response_mat(:,i,:);
end
cell_response_array_norm = cell_response_array./repmat(n_spks(:),[1 length(sfreqCycDeg) length(ori) length(phase0) n_lags]);
%%
freq_response_mat = squeeze(mean(mean(cell_response_array_norm,3),4));
ori_response_mat = squeeze(mean(mean(cell_response_array_norm,2),4));
phase_response_mat = squeeze(mean(mean(cell_response_array_norm,2),3));
emb_ori = [fliplr(-ori(2:end)) ori];
emb_ori_response_mat = cat(2,flipdim(ori_response_mat(:,2:end,:),2),ori_response_mat);


%% FIND RELAVENT TIME LAGS (BEST ARE 1 and 2-sample lags)
freq_mod_index_mat = squeeze((max(freq_response_mat,[],2) - min(freq_response_mat,[],2))./(max(freq_response_mat,[],2) + min(freq_response_mat,[],2)));
%imagesc(freq_mod_index);
[~,best_lag] = max(freq_mod_index_mat,[],2);
for i = 1:n_cells
    [~,max_floc(i)] = max(freq_response_mat(i,:,best_lag(i)),[],2);
    freq_mod_index(i) = freq_mod_index_mat(i,best_lag(i));
end
preferred_freqs = sfreqCycDeg(max_floc);
%% COMPUTE ORIENTATION TUNING
%pull out response at best frequency
for i = 1:n_cells
    overall_orientation_tuning(i,:,:,:) = cell_response_array_norm(i,max_floc(i),:,:,:);
end
ph_inv = squeeze(mean(overall_orientation_tuning,3)); %average over phases
ori_mod_index_mat = squeeze((max(ph_inv,[],2) - min(ph_inv,[],2))./(max(ph_inv,[],2) + min(ph_inv,[],2)));
[~,best_lag] = max(ori_mod_index_mat,[],2);
for i = 1:n_cells
    orientation_tuning(i,:) = ph_inv(i,:,best_lag(i));
    ori_mod_index(i) = ori_mod_index_mat(i,best_lag(i));
end
[~,best_orientation] = max(orientation_tuning,[],2);
preferred_orientation = ori(best_orientation);
%% ALIGN ORIENTATION TUNING
orientation_modulation = (max(orientation_tuning,[],2) - min(orientation_tuning,[],2))./...
    (max(orientation_tuning,[],2) + min(orientation_tuning,[],2));
for i = 1:n_cells
    aligned_orientation_tuning(i,:) = circshift(orientation_tuning(i,:),[0 -best_orientation(i)+9]);
end
for i = 1:n_cells
    phase_dep(i,:) = overall_orientation_tuning(i,best_orientation(i),:,best_lag(i));
end
phase_mod_index = (max(phase_dep,[],2) - min(phase_dep,[],2))./(max(phase_dep,[],2) + min(phase_dep,[],2));
[~,preferred_phase] = max(phase_dep,[],2);
preferred_phase = phase0(preferred_phase);

cd ~/James_scripts/GLM/t1/
save hartley_stats freq_mod_index preferred_freqs ori_mod_index preferred_orientation phase_mod_index preferred_phase