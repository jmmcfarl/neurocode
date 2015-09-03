function [scale,offset,SSE,SSE2,DC] = get_temp_fit(dataseg,template)

scale = (dot(dataseg,template)-sum(template)*sum(dataseg)/length(dataseg))...
    /(sum(template.^2)-sum(template)*sum(template)/length(template));

offset = (sum(dataseg)-scale*sum(template))/length(template);

SSE = sum(dataseg.^2)+scale^2*sum(template.^2)+length(template)*offset^2 ...
    -2*(scale*dot(dataseg,template)+offset*sum(dataseg)-scale*offset*sum(template));
SSE = sqrt(SSE/(length(template)-1));

fitfun = scale*template + offset;
SSE2 = sum((fitfun'-dataseg).^2);
SSE2 = sqrt(SSE2/(length(template)-1));
DC = scale/SSE;