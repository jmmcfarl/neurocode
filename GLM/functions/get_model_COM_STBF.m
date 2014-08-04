function model = get_model_COM_STBF(model,silent)

if nargin < 2
    silent = 1;
end

sdim = model.mods(1).SDIM;
nmods = length(model.mods);
kern_len = length(model.mods(1).k);
kern_t = kern_len/sdim;
STCbvs = model.STCbasis;

x_ax = 1:sdim;

for i = 1:nmods
    spatial_dist = var(reshape(STCbvs*model.mods(i).STCcf,kern_t,sdim));
    spatial_dist = spatial_dist/sum(spatial_dist);
%     [maxval ,maxloc] = max(spatial_dist);
%     model.mods(i).COM = x_ax(maxloc);
    model.mods(i).COM = spatial_dist*x_ax';
end

if ~silent
   for i = 1:nmods
      fprintf('Module %d, COM: %.5f\n',i,model.mods(i).COM);
   end
end