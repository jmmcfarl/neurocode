dt   = 19.9920031987/1000;
% dt   = 20/1000;
sdim = 64;

rstimfiles = {'ara-64p-50h-2m-mn125-ct015.s0', 'wra-64p-50h-2m-mn125-ct015.s1', 'pra-64p-50h-2m-mn125-ct015.s2', 'nra-64p-50h-2m-mn125-ct015.s3', ...
    'nra-64p-50h-2m-mn125-ct025.s4', 'ara-64p-50h-2m-mn125-ct025.s5', 'wra-64p-50h-2m-mn125-ct025.s6', 'pra-64p-50h-2m-mn125-ct025.s7',...
    'pra-64p-50h-2m-mn125-ct035.s8', 'ara-64p-50h-2m-mn125-ct035.s9', 'nra-64p-50h-2m-mn125-ct035.s10', 'wra-64p-50h-2m-mn125-ct035.s11',...
    'nra-64p-50h-2m-mn125-ct045.s12', 'ara-64p-50h-2m-mn125-ct045.s13', 'pra-64p-50h-2m-mn125-ct045.s14', 'wra-64p-50h-2m-mn125-ct045.s15'};

stimfiles  = cellfun(@(x)x(1:26),rstimfiles,'UniformOutput',0);
nstims = length(stimfiles);

datdir = '~/Data/blanche/rec_74';
%stimfiles = dir([datdir,'/stim_data/*.stim']);
cd(datdir)
cd matlabdata/
load ./dstimpsRec74
load ./spksegsRec74
used_filts = [...
    1 0 0 1 0;
    1 0 0 0 0;
    1 1 0 0 1;
    1 1 0 1 0;
    1 0 1 1 0;
    1 1 0 0 0;
    0 1 1 0 0;
    1 0 0 0 0;
    1 0 0 0 0;
    1 0 0 1 0;
    1 0 0 0 0;
    1 1 1 0 0;
    1 1 1 1 0;
    1 0 0 0 0;
    1 0 0 0 1;
    1 0 1 0 0;
    1 0 0 0 0;
    1 0 0 0 0;
    1 1 0 0 0;
    1 1 0 0 0;
    1 0 0 0 0;
    0 1 0 0 0;
    1 0 1 0 0;
    1 0 0 0 0;
    1 0 0 0 0;
    1 1 1 1 1
    1 0 0 0 1];

%%

for mod_cell_num = 1:27;
    mod_cell_num
    cd '/Users/James/James_scripts/GLM/t1/allcell_fits/'
    eval(sprintf('load cell%d_500wdims_5mods_350dX_40L1',mod_cell_num));
%     plot2d_mod(full_glm,0)
    
    used_mods = find(used_filts(mod_cell_num,:)==1);
    full_glm.mods = full_glm.mods(used_mods);
%     plot2d_mod(full_glm,0)
    for i = 1:length(full_glm.mods)
        full_glm.mods(i).nltype = 'uncon';
    end
    full_glm.basis = 'pix';
    %%
    stim_set = 9;
    
    flen = 6;
    S0 = makeStimRows(dstimps74{stim_set},flen);
    S0 = S0/std(S0(:));
    %%
    % rep_stim = repmat(dstimps74{stim_set},24,1);
    % Srep = makeStimRows(rep_stim,flen)/std(S0(:));
    % spkbs = 1+floor(spksegs74{mod_cell_num,stim_set}/dt);
    % pix_mat = get_pix_mat(full_glm);
    [full_glm,norm_vals] = normalizeRFs_full(full_glm,S0);
    % full_glm2 = fitWeights_full(full_glm,Srep*pix_mat,spkbs,1);
    
    %%
    pred_rate = get_predicted_rate(full_glm,S0)/dt;
    t_axis = dt*(1:250);
    %%
    for used_cell_num = 1:25;
        mod_spikes = mod(spksegs74{used_cell_num,stim_set},250*dt);
        psth = hist(mod_spikes,t_axis)/dt;
        [a,b] = corrcoef(psth,pred_rate);
        mod_corr(mod_cell_num,used_cell_num) = a(2,1);
        mod_corr_p(mod_cell_num,used_cell_num) = b(2,1);
%         raster(mod_spikes);
%         hold on
%         plot(t_axis,(pred_rate-mean(pred_rate))/std(pred_rate)*3-20,'r')
%         plot(t_axis,(psth-mean(psth))/std(psth)*3-20,'k')
%         
%         pause
%         close all
    end
end