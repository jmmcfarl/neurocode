function TT_times = load_t_file(tfilelist)

S = LoadSpikes(tfilelist);
num_files = length(S);
for i = 1:length(S);
    TT_times{i} = Data(S{i});
    TT_times{i} = TT_times{i}*100;
end
