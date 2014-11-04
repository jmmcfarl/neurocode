clear all

cd C:\WC_Germany\persistent_downs\
load overall_EC_dir sess_data 
parietal = find(strcmp({sess_data(:).region},'parietal'));
frontal = find(strcmp({sess_data(:).region},'frontal'));
prefrontal = find(strcmp({sess_data(:).region},'prefrontal'));
barrel = find(strcmp({sess_data(:).region},'barrel'));
v1 = find(strcmp({sess_data(:).region},'V1'));

pyramidal = find(strcmp({sess_data(:).cell_type},'pyramidal'));

% L23 = strcmp({sess_data(:).layer},'23');
% L5 = strcmp({sess_data(:).layer},'5');
% L56 = strcmp({sess_data(:).layer},'56');
% L34 = strcmp({sess_data(:).layer},'34');

data_cnt = 0;

for ii = 1:length(parietal)
    data_cnt = data_cnt + 1;
    data(data_cnt).id = data_cnt;
    data(data_cnt).dir = sess_data(parietal(ii)).directory;
    data(data_cnt).region = 'parietal';
    data(data_cnt).cell_type = sess_data(parietal(ii)).cell_type;
    data(data_cnt).layer = sess_data(parietal(ii)).layer;
    data(data_cnt).ep = Inf;
    data(data_cnt).dp = Inf;
end
for ii = 1:length(frontal)
    data_cnt = data_cnt + 1;
    data(data_cnt).id = data_cnt;
    data(data_cnt).dir = sess_data(frontal(ii)).directory;
    data(data_cnt).region = 'frontal';
    data(data_cnt).cell_type = sess_data(frontal(ii)).cell_type;
    data(data_cnt).layer = sess_data(frontal(ii)).layer;
    data(data_cnt).ep = Inf;
    data(data_cnt).dp = Inf;
end
for ii = 1:length(prefrontal)
    data_cnt = data_cnt + 1;
    data(data_cnt).id = data_cnt;
    data(data_cnt).dir = sess_data(prefrontal(ii)).directory;
    data(data_cnt).region = 'prefrontal';
    data(data_cnt).cell_type = sess_data(prefrontal(ii)).cell_type;
    data(data_cnt).layer = sess_data(prefrontal(ii)).layer;
    data(data_cnt).ep = Inf;
    data(data_cnt).dp = Inf;
end
for ii = 1:length(barrel)
    data_cnt = data_cnt + 1;
    data(data_cnt).id = data_cnt;
    data(data_cnt).dir = sess_data(barrel(ii)).directory;
    data(data_cnt).region = 'barrel';
    data(data_cnt).cell_type = sess_data(barrel(ii)).cell_type;
    data(data_cnt).layer = sess_data(barrel(ii)).layer;
    data(data_cnt).ep = Inf;
    data(data_cnt).dp = Inf;
end
for ii = 1:length(v1)
    data_cnt = data_cnt + 1;
    data(data_cnt).id = data_cnt;
    data(data_cnt).dir = sess_data(v1(ii)).directory;
    data(data_cnt).region = 'V1';
    data(data_cnt).cell_type = sess_data(v1(ii)).cell_type;
    data(data_cnt).layer = sess_data(v1(ii)).layer;
    data(data_cnt).ep = Inf;
    data(data_cnt).dp = Inf;
end

%%
cd C:\WC_Germany\final_pdown_analysis\
save compiled_corticalMP_data data
