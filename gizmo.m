function res = gizmo(dat,varargin)

% res = gizmo(dat,...)
%
% Writes dat as a binary file to disk, prepares a script with a call to R
% functions that will compute model estimation in R.
% 
% optional name, value pairs
%   'fbasename' = basename for all files created (default = 'gizmo').
%   'frame'     = structure to write as data.frame for model.
%   'asfactors' = logical vector indicating which fields of frame should be
%                 treated as factors in R (and which ones not).
%   'zscore'    = cell array of strings indicating which fields should be
%                 zscored before writing the frame.
%   'formula'   = formula string of the model.
%   'Rfun'      = R function to use for modeling.
%   'Rargs'     = other arguments to Rfun (other than frame and formula).
%   'doRun'     = boolean whether to run the model in R or not.
%   'rootdir'   = the root directory in which to run and save
%
%

% reading arguments and setting default values
def = struct('fbasename','unnamed_gizmo',...
    'frame',struct(),...
    'formula','~ ',...
    'asfactors',NaN,...
    'zscore',{NaN},...
    'Rfun','gizlm',...
    'doRun',1,...
    'rootdir',pwd);

if numel(varargin) == 1
    cfg = setdef(varargin{1},def);
else
    cfg = setdef(vararg2struct(varargin),def);
end
if iscellstr(cfg.asfactors)
    str = 'c(';
    for i = 1:numel(cfg.asfactors)
        str = [str '"' cfg.asfactors{i} '",'];
    end
    str(end) = ')';
    cfg.asfactors = str;
end

defRargs = struct('Rargs',struct(...
    'formula', cfg.formula, ...
    'basename', ['"' cfg.fbasename '"'],...
    'nblocks', 1000,...
    'asfactors',cfg.asfactors));
cfg = setdef(cfg,defRargs);

cfg.Rargs = struct2vararg(cfg.Rargs);

goto(cfg.rootdir);
% first write binary data to disk
fid = fopen(fullfile(cfg.rootdir,[cfg.fbasename '.dat']),'w','ieee-le');
if fid == -1
    error('Cannot write file. Check permissions and space.')
end
count = fwrite(fid,dat,'single');
fclose(fid);
if count ~= numel(dat)
    error('Error writing data to file');
else wdatok = 1;
end

% then write frame to txt
% first zscore if required
if iscell(cfg.zscore) || not(any(isnan(cfg.zscore)))
    fs = fieldnames(cfg.frame);
    if iscell(cfg.zscore)% assume we're naming fields to zscore.
        toz = cellfun(@(x)any(strcmp(x,cfg.zscore)),fs);
        cfg.zscore = toz;
    end
    for i = 1:numel(fs)
        if cfg.zscore(i)
            if isnumeric(cfg.frame(1).(fs{i}))
                [cfg.frame.(fs{i})] = rep2struct(zscore([cfg.frame.(fs{i})]));
            else
                error('Cannot zscore non numeric predictor');
            end
        end
    end
end

switch class(cfg.frame)
    case 'table'
        writetable(cfg.frame,fullfile(cfg.rootdir,[cfg.fbasename '.df']),'delimiter','\t','filetype','text');
        wfrok = 1;
    case 'struct'
        wfrok = write_table(fullfile(cfg.rootdir,[cfg.fbasename '.df']),mystruct2table(cfg.frame));
    otherwise
        error('unknown type for frame')
end

% write a little script for R
fid = fopen(fullfile(cfg.rootdir,[cfg.fbasename '.R']),'wt');
rfuns = which('gizfuns.R');
if isempty(rfuns)
    error('could not find R functions')
else
    str = ['source("' rfuns '")'];% that loads giz R functions
end
fprintf(fid,'%s\n',str);
str = ['setwd("' cfg.rootdir '")'];
fprintf(fid,'%s\n',str);

str = [cfg.Rfun '( '];
for i = 1:2:numel(cfg.Rargs)
    if ischar(cfg.Rargs{i+1})
        str = [str cfg.Rargs{i} ' = ' cfg.Rargs{i+1}];
    elseif isnumeric(cfg.Rargs{i+1}) || islogical(cfg.Rargs{i+1})
        str = [str cfg.Rargs{i} ' = ' num2c(cfg.Rargs{i+1})];
    elseif iscellstr(cfg.Rargs{i+1})
        str = [str cfg.Rargs{i} ' = ' cstr2c(cfg.Rargs{i+1})];
    else
        str = [str cfg.Rargs{i} ' = ' num2str(cfg.Rargs{i+1})];
    end
    if i+1 == numel(cfg.Rargs)
        str(end+1) = ')';
    else
        str(end+1:end+2) =  ', ';
    end
end
fprintf(fid,'%s\n',str);
fclose(fid);

% write a batch that will call R and run the above script
fid = fopen('Runscript.bat','wt');
str = ['R  CMD BATCH --no-save --no-restore ' [cfg.fbasename '.R'] ];
fprintf(fid,'%s\n',str);
fclose(fid);

% delete any previously produced data generated with this fbasename
delete([cfg.fbasename '_*.dat'])
delete([cfg.fbasename '_*.txt'])
% create info.txt
fid = fopen(fullfile(cfg.rootdir,[cfg.fbasename '_info.txt']),'wt');
fprintf(fid,'data size: ');
fprintf(fid,'%d\t',size(dat));
fprintf(fid,'\n');
fclose(fid);

if isunix
    !chmod +x Runscript.bat
end

if cfg.doRun
    !./Runscript.bat
    
    res = giz_readfiles(cfg.fbasename,size(dat));
else
    res = 1;
end

goback;

return


function m = mat2cells(c)

str = 'mat2cell(c,';
for i = 1:ndims(c)
    str = [str 'ones(size(c,' num2str(i) '),1),'];
end
str(end) = ')';
m = eval(str);

function out = num2c(in)
% create a string c(in(1),in(2)...) to pass a vector to R
out = 'c(';
for i = 1:numel(in)
    out = [out num2str(in(i)) ','];
end
out(end) = ')';

function out = cstr2c(in)
% create a string c(in(1),in(2)...) to pass a vector to R
out = 'c(';
for i = 1:numel(in)
    out = [out '"' in{i} '",'];
end
out(end) = ')';


function txt = mystruct2table(T,header,fields,sname)

% txt = struct2table(T,header,fields,sname)
% Turn fields fields of structure T into a cell array "txt".
% If header is true, a first header line is added with field names. If
% header is a cell array of strings, these strings are used as headers
% (length(header) must be == length(fields))
% If provided, sname is a string that will be added as a first column to
% each row of txt.

% if nargin == 0
%     T = evalin('caller','T');
%     fields = {'RefDur' 'Freq' 'RT' 'RealT' 'FB'};
% end
if not(exist('fields','var'))
    fields = fieldnames(T);
end
if not(exist('header','var'))
    header = true;
end
if iscellstr(header)
    headerline = header;
    header = true;
else
    headerline = fields;
end
if not(numel(headerline) == numel(fields))
    error('Number of columns inconsistency, check input');
end
gotsname = exist('sname','var');
txt = cell(numel(T)+header,numel(fields)+gotsname);
if header
    if gotsname
        txt(1,:) = {'suj' headerline{:}};
    else
        txt(1,:) = {headerline{:}};
    end
end

for itri = 1:numel(T)
    if gotsname
        txt{itri+header,1} = sname;
    end
    for i_f = 1:numel(fields)
        txt{itri+header,i_f+gotsname} = T(itri).(fields{i_f});
    end
end







