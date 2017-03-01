
addpath(cdhome('nuGiZmo'))

opts.rmchan = {'EXG.*'};
opts.timewindow = [-200 1000]/1000;
opts.reref = 'EXG[45]';
opts.nboot = 500;


EEG = pop_loadset('/Osz_04/Max/N399Ctxt/data/S008_strength/S008_strength_viscor_icacor.set');
if isfield(opts,'reref')
    EEG = pop_reref(EEG,chnb(opts.reref),'keepref','on');
end
if isfield(opts,'rmchan')
    EEG = pop_select( EEG,'nochannel',chnb(opts.rmchan));
end
if isfield(opts,'timewindow')
    EEG = pop_select(EEG,'time',opts.timewindow);
end

%%
tomodel = permute(EEG.data,[3 2 1]);% [trials,time, channels]

EOI = 130;
FOI = {'RelatednessCont','consensus_partner','consensus_prime'};
event = EEG.event;
event = event(strcmp({event.type},num2str(EOI)));
fs = fieldnames(event);
event = rmfield(event,setxor(fs,FOI));

results = gizmo(tomodel,'fbasename','S008nunutry'...
    ,'frame',event...
    ,'asfactors',false(1,numel(fieldnames(event)))...
    ,'zscore',true(1,numel(fieldnames(event)))...
    ,'formula','~ RelatednessCont + consensus_partner + consensus_prime + RelatednessCont:consensus_partner + RelatednessCont:consensus_prime + consensus_partner:consensus_prime'...
    ,'doRun',1 ...
    ,'Rfun','gizfastlm');
%% if not doRun
% tomodel = giz_loadbin('S008firstry.dat');
% tomodel = reshape(tomodel,[372,1024,69]);

results2 = giz_readfiles('S008nunutry',size(tomodel));
figure;imagesc(results.design_dat)

figure;
[tpts, tvals]= timepts([-200 1000]);
imagesc(tvals,1:EEG.nbchan,permute(results.coefs_dat(1,tpts,:),[3 2 1]));
set(gca,'clim',[-15 15])

%%
s = size(tomodel(:,:));
b = zeros(size(results.design_dat,2),s(2),'single');
h = waitbar(0);
pp = pinv(results.design_dat);
for i = 1:s(2)
    if ~rem(i,100)
        waitbar(i/s(2),h,num2str(i));
    end
    b(:,i) = pp * tomodel(:,i);
end
close(h)
%%
%all(b(:) == results.coefs_dat(:))
b = reshape(b,size(results.coefs_dat));
figure; 
[tpts tvals]= timepts([-200 1000]);
imagesc(tvals,1:EEG.nbchan,permute(b(1,tpts,:),[3 2 1]));
set(gca,'clim',[-15 15])

