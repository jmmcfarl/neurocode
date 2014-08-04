function [template_scores,templates_used] = get_template_scores(X,templates,template_params)

[N_spks,D,n_chs] = size(X);
n_templates = size(templates,2);
channels = template_params.channels;
t_derivs = template_params.t_derivs;
use_ms = template_params.mean_sub;
if nargin < 3 || isempty(channels)
    channels = ones(n_templates,1);
end
if nargin < 4 || isempty(t_derivs)
   t_derivs = zeros(n_templates,1); 
end
if nargin < 5 || isempty(use_ms)
    use_ms = ones(n_templates,1);
end
if any(t_derivs > 1)
    error('Only first derivative supported');
end

template_scores = nan(N_spks,n_templates);
for tt = 1:n_templates
   cur_spikes = squeeze(X(:,:,channels(tt)));
   if t_derivs(tt) == 1
       cur_spikes = [zeros(N_spks,1) diff(cur_spikes,1,2)];
   end
   if use_ms(tt) == 1
      cur_spikes = bsxfun(@minus,cur_spikes,mean(cur_spikes,2));
      templates(:,tt) = templates(:,tt) - mean(templates(:,tt));
   end
   template_scores(:,tt) = cur_spikes*templates(:,tt);
end
template_scores = bsxfun(@rdivide,template_scores,std(template_scores));