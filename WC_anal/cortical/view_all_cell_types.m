function [] = view_all_cell_types(data_mat, xvals)

cd C:\WC_Germany\Cortical_analysis

load C:\WC_Germany\Cortical_analysis\cortical_dir

dividerlines = [11.5 14.5 24.5 28.5 31.5 36.5 44.5 49.5 54.5 55.5]

minx = min(xvals);
maxx = max(xvals);
imagesc(xvals,1:size(data_mat,1),data_mat);
shading flat
for i = 1:length(dividerlines)
    line([minx maxx],[dividerlines(i) dividerlines(i)],'Color','w')
end