function [res, meta] = giz_readfiles(fbasename,expectedsize)

if not(exist('expectedsize','var')) || isempty(expectedsize)
    ff = rdir([fbasename '_info.txt']);
    if isempty(ff)
        error(['Info file ' fbasename '_info.txt' ' missing']);
    end
    txt = readtext(ff.name);
    fs = regexp(txt,'([^:]+):(.*)','tokens');
    for i_f = 1:numel(fs)
        f = strrep(fs{i_f}{1}{1},' ','_');
        c = str2num(fs{i_f}{1}{2});
        if isempty(c)
            c = fs{i_f}{1}{2};
        end
        meta.(f) = c;
    end
    expectedsize = meta.data_size;
else
    meta.data_size = expectedsize;
end
    

res = struct();
fs = rdir([fbasename '_*.dat']);
for i = 1:numel(fs)
    field = regexp(fs(i).name,[fbasename '_(.*).dat'],'tokens');
    field = [field{1}{1} '_dat'];
    d = loadbin(fs(i).name);
    if numel(d) == prod(expectedsize) % residuals
        res.(field) = reshape(d,expectedsize);
    elseif ~rem(numel(d),prod(expectedsize(2:end))) % coefficients
        s = mat2cells(expectedsize(2:end));
        res.(field) = reshape(d,[],s{:});
    elseif ~rem(numel(d),expectedsize(1)) % design matrix
        res.(field) = reshape(d,expectedsize(1),[]);
    end
end
fs = rdir([fbasename '_*.txt']);
for i = 1:numel(fs)
    field = regexp(fs(i).name,[fbasename '_(.*).txt'],'tokens');
    field = [field{1}{1} '_txt'];
    if strcmp(field,'info_txt')
        continue
    end
    d = readtext(fs(i).name,'\t');
    d = regexprep(d,'["\(\)]','');
    res.(field) = d;
end


function d = loadbin(fn)

fid = fopen(fn,'rb','l');
d = fread(fid,Inf,'single=>single');
fclose(fid);
