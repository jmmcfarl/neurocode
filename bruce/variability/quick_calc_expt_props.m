% clear all
% close all

addpath('~/other_code/fastBSpline/');

global Expt_name bar_ori monk_name rec_type rec_number Elist_cnt bori_cnt

% Expt_name = 'M012';
% monk_name = 'jbe';
% bar_ori = 0; %bar orientation to use (only for UA recs)
% rec_number = 1;
% %
% [266-80 270-60 275-135 277-70 281-140 287-90 289-160 294-40 296-45 297-0/90 5-50 9-0 10-60 11-160 12-0 13-100 14-40 320-100]

%%
data_dir = ['~/Data/bruce/' Expt_name];
if ~exist(data_dir,'dir')
    data_dir = ['/media/NTlab_data3/Data/bruce/' Expt_name];
end
if ~exist(data_dir,'dir');
    error('Couldnt find data directory');
end
anal_dir = ['~/Analysis/bruce/' Expt_name '/variability/'];
if ~exist(anal_dir,'dir')
    mkdir(anal_dir)
end

Edata_file = strcat(data_dir,'/',monk_name,Expt_name,'Expts');
load(Edata_file);

%is this a laminar probe or utah array rec?
ff = find(cellfun(@(x) ~isempty(x),Expts),1);
if strcmp(Expts{ff}.Header.DataType,'GridData 96')
    rec_type = 'UA';
elseif strcmp(Expts{ff}.Header.DataType,'Spike2')
    rec_type = 'LP';
end

%load in packaged data
data_name = sprintf('%s/packaged_data_ori%d',data_dir,bar_ori);
if rec_number > 1
    data_name = strcat(data_name,sprintf('_r%d',rec_number));
end
fprintf('Loading %s\n',data_name);
load(data_name);

%%

out_data(Elist_cnt,bori_cnt).wi_mode = mode(expt_data.expt_wi)/params.scale_fac;
out_data(Elist_cnt,bori_cnt).dw_mode = mode(expt_data.expt_dw)/params.scale_fac;
out_data(Elist_cnt,bori_cnt).fw_mode = mode(cellfun(@(x) x.Stimvals.fw,Expts(expt_data.used_blocks)));
out_data(Elist_cnt,bori_cnt).dds_mode = mode(expt_data.expt_dds);