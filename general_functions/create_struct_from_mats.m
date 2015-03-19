function S = create_struct_from_mats(varargin)

if length(varargin) < 2
    error('Need at least a pair of inputs');
end

j = 1;
field_cnt = 0;
while j <= length(varargin)
field_cnt = field_cnt + 1;

vname{field_cnt} = varargin{j};
vmat = varargin{j+1};
[N,K] = size(vmat);

cellSize = ones(N,K);
vcell{field_cnt} = mat2cell(vmat,cellSize);

j = j + 2;
end

command = sprintf('S = struct(');
for ii = 1:field_cnt
   command = strcat(command,sprintf('''%s'',vcell{%d},',vname{ii},ii)); 
end
command([end end+1]) = ');';
eval(command);