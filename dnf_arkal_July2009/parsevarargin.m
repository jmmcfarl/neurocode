function varstruct = parsevarargin(v, templatestruct)
% varstruct = parsevarargin(v, templatestruct)
%
% INPUTS:
%  v = varargin
%  templatestruct = structure with fieldnames and default values
%
%
% NOTES:
%  fieldnames and argument names are case insensitive
%
% David Nguyen <dpnguyen@mit.edu>
%
% $id$
%
  
  if ~exist('templatestruct')
    [dummy, vars] = parseparms(v);
    fnames = vars(1:2:end);
    vals = vars(2:2:end);
    varstruct = [];
    for k = 1:length(fnames)
      varstruct = setfield(varstruct, fnames{k}, vals{k});
    end
  else
    varstruct = templatestruct(1);
    F = fieldnames(templatestruct);
    f = lower(F);
    
    J = 1;
    while J <= (length(v)-1)
      
      if ischar(v{J})        
        id  = find(strcmp(f, lower(v(J))));
        if length(id) > 0
          varstruct = setfield(varstruct, F{id(1)}, v{J+1});
          J = J + 1;
        end
      end
      J = J + 1;

    end
  end
    
    
    
