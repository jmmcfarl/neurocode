function effect_size = get_effect_size(u,sig)

pooled_s = sqrt(mean([sig(1)^2 sig(2)^2]));

effect_size = abs(u(1)-u(2))/pooled_s;