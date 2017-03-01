addpath(cdhome('nuGiZmo'))

%% one dataset

ntri = 100; ntime = 1000;

% effects = [cond=red dist cond=red:dist];
effects = [2 3 1];

[MEEG,frame] = onesuj(ntri,ntime,effects);

figure(499);clf;
imagesc(MEEG)
%% with gizmo
tic

res = gizmo(MEEG,'frame',frame,'asfactor',{'cond'},'formula', '~ cond * dist');

% this should recover effects (+ intercept)
mean(res.coefs_dat,2)
toc

%% with fitlm (slowest)
tic
tblframe = struct2table(frame);

for i = 1:size(MEEG,2)
    tblframe.Y = MEEG(:,i);
    m = fitlm(tblframe,'Y~ cond * dist');
    res.coefs_dat_fitlm(:,i) = m.Coefficients.Estimate;
end
mean(res.coefs_dat_fitlm,2)
toc
%% with mvregress (faster)
tic
X = squeeze(struct2cell(frame));
X(1,:) = cellfun(@(X)1+(strcmp(X,'blue')),X(1,:),'uniformoutput',0);
X = cell2mat(X)';
X = x2fx(X,'interaction',1);
for i = 1:size(MEEG,2)
    y = MEEG(:,i);
    b = mvregress(X,y);
    res.coefs_dat_mvregress(:,i) = b;
end
mean(res.coefs_dat_mvregress,2)
toc
%% with pinv (even faster)
tic
X = squeeze(struct2cell(frame));
X(1,:) = cellfun(@(X)1+(strcmp(X,'blue')),X(1,:),'uniformoutput',0);
X = cell2mat(X)';
X = x2fx(X,'interaction',1);
pX = pinv(X);
for i = 1:size(MEEG,2)
    y = MEEG(:,i);
    b = pX * y;
    resids = y - (X * b);
    res.coefs_dat_mvregress(:,i) = b;
end
mean(res.coefs_dat_mvregress,2)
toc
%% faster with gizmo

tic

res = gizmo(MEEG,'frame',frame,'asfactor',{'cond'},'formula', '~ cond * dist','Rfun','gizfastlm');

% this should recover effects (+ intercept)
mean(res.coefs_dat,2)
toc


%%


% below I would like to extend to several subjects...
ntri = 100; ntime = 1000;
effects = [2 3 1];

nsuj = 12;
tic
MEEG = {};frame = {}; res = {};
% level 1
% run each subject
for isuj = 1:nsuj
    [MEEG{isuj}, frame{isuj}] = onesuj(ntri,ntime,effects);
    res{isuj} = gizmo(MEEG{isuj},'frame',frame{isuj},'asfactor',{'cond'},'formula', '~ cond * dist','Rfun','gizfastlm');
end
toc
% level 2
% collect data into one matrix
MEEG2 = cellfun(@(x)x.coefs_dat,res,'uniformoutput',0);
MEEG2 = cat(1,MEEG2{:});
%
% create a frame with conditions
frame2 = cellfun(@(x) x.coefnames_txt,res,'uniformoutput',0);
frame2 = cell2table(cat(1,frame2{:}),'variablenames',{'condeffect'});

res2 = gizmo(MEEG2,'frame',frame2,'asfactor','cond','formula', '~ condeffect - 1');
mean(res2.coefs_dat,2)
res2.coefnames_txt