% Delta, theta, traces ofr bp linear, depths and strips separate
% grids, vs strips vs depths 

% what gamma does across grids strips depths 
% do individual lines x axis distance, power y axis
%

% permutation tes twould be on distance/frequency for perumutation windows 
% speech vs. non speech IFG, stim vs non stim STG - make sure its
% exclusive, non speech should not have any stim, non stim should not have
% any speech 

% spike analysis for 4 vs 8 vs 12 
% skip 2/3 rows 



%bipolar vs referntial spectrum
% take chanel 1 and 2 , take mean of individual spectrum 
% take bipolrar reference 1 vs 2 , and then spectrum, and compariuson of
% this adnt he above 
% then plot the difference between the two lines (bipolar specturm vs
% linear) 

% better represetnation of what youre doing vs referential 


% microgrid 
% < 10 mm aspect 



% BIPOLAR PAIR ANALYSIS: EACH VS. ALL

savePlots = true; 
folderFigures = '/Users/davidcaldwell/Box/Kleenlab/David/Results/Apr Results/';
pt='EC175'; % EC175 and EC183 both have intact 16x16 square grids (channel #s 1:256)
nchtocheck=256; %***ideally, should eventually convert this into an index of channels to check (STG, IFG, etc) which could be defined before the spectra loop to save computational time...
windowstocheck=1:100;
binsz=3; % bin size in mm

onlygrids=true;
onlydepths=false;
onlystrips=false;

cm=cool(6); cm(1,:)=[0 0 0];
datadir='/Volumes/KLEEN_DRIVE/David/Bipolar project/baseline-high-density-data/'; %bandpassfiltered/';
cd([datadir])
load('/Volumes/KLEEN_DRIVE/David/Bipolar project/taggedspikes.mat')
sfx=512;
frxrange=[2 200]; %frequency range to examine
ft=[2 5 10 20 50 100 200]; ftl=cellstr(num2str(ft')); %frequency labels for plots

u=dir; uptbl={}; for i=1:length(u); uname=u(i).name; uptbl{i,1}=uname(1:end-28); end; uptbl(1:2)=[]; clear i u uname

p=find(strcmpi(pts,pt)); %patient number ID
pblocks=strfind(uptbl,pts{p});
for i=1:length(pblocks);
    isbl(i,1)=~isempty(pblocks{i});
end
ptbl=find(isbl); if ~isempty(ptbl); disp(['Loading ' pts{p} ' blocks...']); end

%make vector of stimuli/speech
hasStimTotal = [];
hasSpeechTotal = [];

% load all blocks for this patient and stack their baseline windows together
d=[]; nwind=0;
for b=1:length(ptbl); disp(uptbl{ptbl(b)})
    % load using using "_jk" versions of baseline windows, updated 2/2022
    load([datadir uptbl{ptbl(b)} '_baselineWindows_fromraw.mat'])
    % get rid of baseline windows containing spikes or artifact
    spksarti=hasspk | hasarti;
    nonspks_windows(:,spksarti)=[];
    hasstim(spksarti)=[]; %update indices for which windows overlap with stimuli/speech
    hasspeech(spksarti)=[];
    hasStimTotal = [hasStimTotal hasstim];
    hasSpeechTotal = [hasSpeechTotal hasspeech];
    
    clear hasspkvec hasspk hasartivec hasarti spksarti % now clear spike- and artifact-related variables from workspace
    
    % convert to 3D matrix, combine all windows from consecutive blocks for each patient
    for i=1:size(nonspks_windows,2);
        d(:,:,i+nwind)=nonspks_windows{2,i}';
    end
    nwind=size(d,3);
    
    clear nonspks_windows info
end; clear b

hasStimTotal = logical(hasStimTotal);
hasSpeechTotal = logical(hasSpeechTotal);

nch=size(d,2);

% load electrode component infor (grid/strip/depth and how many linear contacts they have in a row
[bpN,bpT]=xlsread(['/Volumes/KLEEN_DRIVE/David/Bipolar project/AN_ElectrodeInfoTDT.xlsx'],pts{p});
[em,eleclabels,anatomy]=getelecs(pts{p},2);

cm=cool(6); cm(1,:)=[0 0 0];
%datadir='/Volumes/KLEEN_DRIVE/David/Bipolar project/baseline-high-density-data/bandpassfiltered/';
%cd([datadir])
load('/Volumes/KLEEN_DRIVE/David/Bipolar project/taggedspikes.mat')
sfx=512;
frxrange=[2 200]; %frequency range to examine
ft=[2 5 10 20 50 100 200]; ftl=cellstr(num2str(ft')); %frequency labels for plots


% if wanting to only look at grids, strips, or depths, then nan the others
if onlygrids||onlystrips||onlydepths;
    for r=1:size(bpT,1)
        if any(strcmpi(bpT(r,2),{'grid','minigrid'})) && ~onlygrids  || ...
                strcmpi(bpT(r,2),'strip')              && ~onlystrips || ...
                strcmpi(bpT(r,2),'depth')              && ~onlydepths;
            d(:,bpN(r,1):bpN(r,2),:)=nan;
        end
    end; clear r
end


% bad channels
badchI=isnan(mean(mean(d,1),3))'; %any channel with nans in any window will be a "bad channel"
badchI(badchidx{p})=1; %follow up with all marked as bad channels in original preprocessing
x=[]; xbch=true(size(badchI)); for i=1:size(bpT,1); x=[x bpN(i,1):bpN(i,2)]; end; xbch(x)=false;
badchI(xbch)=1; %find and nan empty channels unaccounted for by component rows
d(:,badchI,:)=nan;
okc=~badchI; clear x xbch


%*at this point, isolate speech or stim or non-speech/stim

%*at this point, isolate STG or IFG channels
[STG]=getelecs_region(pt,'stg',2);
[IFG]=getelecs_region(pt,{'po','pt'},2);


%% ALL PAIRS (each vs. all others) analysis and example plot
% ***hint hint: opportunity here to select speech or stim windows!

clear Straces_allch;




dSpeech = d(:,:,(hasSpeechTotal & ~hasStimTotal));
dNoSpeechNoStim = d(:,:,(~hasSpeechTotal & ~hasStimTotal));
dStim = d(:,:,(hasStimTotal & ~hasSpeechTotal));

dSpeech = dSpeech(:,:,windowstocheck);
dNoSpeechNoStim = dNoSpeechNoStim(:,:,windowstocheck);
dStim = dStim(:,:,windowstocheck);

[mSpeech,Mbpdist,frx]=bpspectra_EachVsAll(dSpeech,sfx,frxrange,em,nchtocheck);
[mNoSpeechNoStim,Mbpdist,frx]=bpspectra_EachVsAll(dNoSpeechNoStim,sfx,frxrange,em,nchtocheck);
[mStim,Mbpdist,frx]=bpspectra_EachVsAll(dStim,sfx,frxrange,em,nchtocheck);

%mSpeech = M(:,:,:,(hasSpeechTotal & ~hasStimTotal));
%mNoSpeechNoStim = M(:,:,:,(~hasSpeechTotal & ~hasStimTotal));
%mStim = M(:,:,:,(hasStimTotal & ~hasSpeechTotal));

d=d(:,:,windowstocheck); clear Straces_allch; %free up RAM by getting rid of whatever won't be used (only using first ___ number of windows)
[M,Mbpdist,frx]=bpspectra_EachVsAll(d,sfx,frxrange,em,nchtocheck);




%% mean across windows
M=sq(mean(M,4));

%% now that we hav ebipolar pairs and spectra, we can select the channels we are interested in 

% which channnels are we interseted in
chansIntSTG = STG;
chansIntSTG = chansIntSTG(chansIntSTG<=256);

mSpeechSTG = mSpeech(chansIntSTG,chansIntSTG,:,:);
mStimSTG = mStim(chansIntSTG,chansIntSTG,:,:);
mNoSpeechNoStimSTG = mNoSpeechNoStim(chansIntSTG,chansIntSTG,:,:);
MbpdistSTG = Mbpdist(chansIntSTG, chansIntSTG);

chansIntIFG = IFG;
chansIntIFG = chansIntIFG(chansIntIFG<=256);

mSpeechIFG = mSpeech(chansIntIFG,chansIntIFG,:,:);
mStimIFG = mStim(chansIntIFG,chansIntIFG,:,:);
mNoSpeechNoStimIFG = mNoSpeechNoStim(chansIntIFG,chansIntIFG,:,:);
MbpdistIFG = Mbpdist(chansIntIFG,chansIntIFG);


%% unstack and line up into 2D matrix to create channels^2 X frequencies for easier indexing-->binning
Mflat=[];
Mflat_bpdist=[];
for c2=1:nchtocheck;
    Mflat=[Mflat; sq(M(:,c2,:))];
    Mflat_bpdist=[Mflat_bpdist; Mbpdist(:,c2)]; % corresponding distance index
end

binz=0:binsz:85; % can change binsz PRN at this point to test different bin resolutions on the plots below
clear mx my mz
nbinz=length(binz)-1;
mz=[];
%bin by distance and take the mean, creating: frequency X distance (binned)
for i=1:nbinz
    mz(:,i)=nanmean(Mflat(Mflat_bpdist>=binz(i) & Mflat_bpdist<binz(i+1),:),1);
    binindex_min_ltnmax(i,:)=[binz(i) binz(i+1)];
end

% Notes:
% mz is frequency X distance (binned), after having averaged across windows
% binindex_min_ltmax tells you for each column of mz what was the
%         1] minimum (>0) distance, and
%         2] the "less than max" (<) distance
%         that were used to index bipolar pairs for that bin
%         Note: first bin includes zero, which corresponds to
%         the bin containing the referential channels

% zscore the log transformed power according to frequency
nfrx=length(frx);
for i=1:nfrx; nns=~isnan(mz(i,:));
    mz_z(i,nns)=zscore(log(mz(i,nns)));
end;

%% unstack and line up into 3D matrix to create channels^2 X frequencies x trials for easier indexing-->binning
[mz_zSpeech,mzSpeech,mflatSpeech,bpdistSpeech ] = bin_zscore_trial(mSpeech,nchtocheck,Mbpdist,binsz,frx);
[mz_zStim,mzStim,mflatStim,bpdistStim ] = bin_zscore_trial(mStim,nchtocheck,Mbpdist,binsz,frx);
[mz_zNoST,mzNoST,mflatNoST,bpdistNoST ] = bin_zscore_trial(mNoSpeechNoStim,nchtocheck,Mbpdist,binsz,frx);


[mz_zSpeechSTG,mzSpeechSTG,mflatSpeechSTG,bpdistSpeechSTG ] = bin_zscore_trial(mSpeechSTG,length(chansIntSTG),MbpdistSTG,binsz,frx);
[mz_zStimSTG,mzStimSTG,mflatStimSTG,bpdistStimSTG ] = bin_zscore_trial(mStimSTG,length(chansIntSTG),MbpdistSTG,binsz,frx);
[mz_zNoSTSTG,mzNoSTSTG,mflatNoSTSTG,bpdistNoSTSTG ] = bin_zscore_trial(mNoSpeechNoStimSTG,length(chansIntSTG),MbpdistSTG,binsz,frx);

[mz_zSpeechIFG,mzSpeechIFG,mflatSpeechIFG,bpdistSpeechIFG ] = bin_zscore_trial(mSpeechIFG,length(chansIntIFG),MbpdistIFG,binsz,frx);
[mz_zStimIFG,mzStimIFG,mflatStimIFG,bpdistStimIFG ] = bin_zscore_trial(mStimIFG,length(chansIntIFG),MbpdistIFG,binsz,frx);
[mz_zNoSTIFG,mzNoSTIFG,mflatNoSTIFG,bpdistNoSTIFG ] = bin_zscore_trial(mNoSpeechNoStimIFG,length(chansIntIFG),MbpdistIFG,binsz,frx);

%%
figure
subplot(1,3,1)
pcolor(binz(2:size(mz_zNoST,3)+1),frx,squeeze(nanmean(mz_zNoST,2))); shading flat; set(gca,'ydir','normal'); ylabel('Frequency (Hz)'); xlabel('Distance (mm)'); set(gca,'fontsize',14); colorbar;
text(max(xlim)+diff(xlim)/4,mean(ylim),'ln(power)','fontsize',12,'rotation',90,'horizontalalignment','center')
title('No Stimulus/Speech (z-scored by frequency)','fontweight','normal')
set(gca,'yscale','log','ytick',ft,'yticklabel',ftl);

subplot(1,3,2)
pcolor(binz(2:size(mz_zSpeech,3)+1),frx,squeeze(nanmean(mz_zSpeech,2))); shading flat; set(gca,'ydir','normal'); ylabel('Frequency (Hz)'); xlabel('Distance (mm)'); set(gca,'fontsize',14); colorbar;
text(max(xlim)+diff(xlim)/4,mean(ylim),'ln(power)','fontsize',12,'rotation',90,'horizontalalignment','center')
title('Speech (z-scored by frequency)','fontweight','normal')
set(gca,'yscale','log','ytick',ft,'yticklabel',ftl);

subplot(1,3,3)
pcolor(binz(2:size(mz_zStim,3)+1),frx,squeeze(nanmean(mz_zStim,2))); shading flat; set(gca,'ydir','normal'); ylabel('Frequency (Hz)'); xlabel('Distance (mm)'); set(gca,'fontsize',14); colorbar;
text(max(xlim)+diff(xlim)/4,mean(ylim),'ln(power)','fontsize',12,'rotation',90,'horizontalalignment','center')
title('Stimulus (z-scored by frequency)','fontweight','normal')
set(gca,'yscale','log','ytick',ft,'yticklabel',ftl);

%%  permutation testing


%for chan = 1:size(mSpeech,1)
[sizeX,~,sizeY] = size(mz_zSpeech);
[clustersSpeech, pValuesSpeech, tSumsSpeech, permutationDistributionSpeech] = permutest(permute(mz_zSpeech,[1,3,2]),permute(mz_zNoST,[1,3,2]),false,[],[],1);
[clustersStim, pValuesStim, tSumsStim, permutationDistributionStim] = permutest(permute(mz_zStim,[1,3,2]),permute(mz_zNoST,[1,3,2]),false,[],[],1);
[clustersStimSpeech, pValuesStimSpeech, tSumsStimSpeech, permutationDistributionStimSpeech] = permutest(permute(mz_zStim,[1,3,2]),permute(mz_zSpeech,[1,3,2]),false,[],[],1);

%for chan = 1:size(mSpeech,1)
[sizeXSTG,~,sizeYSTG] = size(mz_zSpeechSTG);
[clustersSpeechSTG, pValuesSpeechSTG, tSumsSpeechSTG, permutationDistributionSpeechSTG] = permutest(permute(mz_zSpeechSTG,[1,3,2]),permute(mz_zNoSTSTG,[1,3,2]),false,[],[],1);
[clustersStimSTG, pValuesStimSTG, tSumsStimSTG, permutationDistributionStimSTG] = permutest(permute(mz_zStimSTG,[1,3,2]),permute(mz_zNoSTSTG,[1,3,2]),false,[],[],1);
[clustersStimSpeechSTG, pValuesStimSpeechSTG, tSumsStimSpeechSTG, permutationDistributionStimSpeechSTG] = permutest(permute(mz_zStimSTG,[1,3,2]),permute(mz_zSpeechSTG,[1,3,2]),false,[],[],1);

%for chan = 1:size(mSpeech,1)
[sizeXIFG,~,sizeYIFG] = size(mz_zSpeechIFG);
[clustersSpeechIFG, pValuesSpeechIFG, tSumsSpeechIFG, permutationDistributionSpeechIFG] = permutest(permute(mz_zSpeechIFG,[1,3,2]),permute(mz_zNoSTIFG,[1,3,2]),false,[],[],1);
[clustersStimIFG, pValuesStimIFG, tSumsStimIFG, permutationDistributionStimIFG] = permutest(permute(mz_zStimIFG,[1,3,2]),permute(mz_zNoSTIFG,[1,3,2]),false,[],[],1);
[clustersStimSpeechIFG, pValuesStimSpeechIFG, tSumsStimSpeechIFG, permutationDistributionStimSpeechIFG] = permutest(permute(mz_zStimIFG,[1,3,2]),permute(mz_zSpeechIFG,[1,3,2]),false,[],[],1);


%%

[clusterSigSpeech,pValsSigSpeech,boundarySigSpeech,clustXSpeech,clustYSpeech] = signif_boundary(clustersSpeech,pValuesSpeech,sizeX,sizeY);
[clusterSigStim,pValsSigStim,boundarySigStim,clustXStim,clustYStim] = signif_boundary(clustersStim,pValuesStim,sizeX,sizeY);
[clusterSigStimSpeech,pValsSigStimSpeech,boundarySigStimSpeech,clustXStimSpeech,clustYStimSpeech] = signif_boundary(clustersStimSpeech,pValuesStimSpeech,sizeX,sizeY);

[clusterSigSpeechSTG,pValsSigSpeechSTG,boundarySigSpeechSTG,clustXSpeechSTG,clustYSpeechSTG] = signif_boundary(clustersSpeechSTG,pValuesSpeechSTG,sizeXSTG,sizeYSTG);
[clusterSigStimSTG,pValsSigStimSTG,boundarySigStimSTG,clustXStimSTG,clustYStimSTG] = signif_boundary(clustersStimSTG,pValuesStimSTG,sizeXSTG,sizeYSTG);
[clusterSigStimSpeechSTG,pValsSigStimSpeechSTG,boundarySigStimSpeechSTG,clustXStimSpeechSTG,clustYStimSpeechSTG] = signif_boundary(clustersStimSpeechSTG,pValuesStimSpeechSTG,sizeXSTG,sizeYSTG);

[clusterSigSpeechIFG,pValsSigSpeechIFG,boundarySigSpeechIFG,clustXSpeechIFG,clustYSpeechIFG] = signif_boundary(clustersSpeechIFG,pValuesSpeechIFG,sizeXIFG,sizeYIFG);
[clusterSigStimIFG,pValsSigStimIFG,boundarySigStimIFG,clustXStimIFG,clustYStimIFG] = signif_boundary(clustersStimIFG,pValuesStimIFG,sizeXIFG,sizeYIFG);
[clusterSigStimSpeechIFG,pValsSigStimSpeechIFG,boundarySigStimSpeechIFG,clustXStimSpeechIFG,clustYStimSpeechIFG] = signif_boundary(clustersStimSpeechIFG,pValuesStimSpeechIFG,sizeXIFG,sizeYIFG);

%%
avgNoST = nanmean(mz_zNoST,2);
avgSpeech = nanmean(mz_zSpeech,2);
avgStim = nanmean(mz_zStim,2);
maxZ = max([avgNoST(:);avgSpeech(:);avgStim(:)]);
minZ = min([avgNoST(:);avgSpeech(:);avgStim(:)]);

figure
subplot(3,1,1)
pcolor(binz(2:size(mz_zNoST,3)+1),frx,squeeze(nanmean(mz_zNoST,2))); shading flat; set(gca,'ydir','normal'); ylabel('Frequency (Hz)'); xlabel('Distance (mm)'); set(gca,'fontsize',14); colorbar;
text(max(xlim)+diff(xlim)/4,mean(ylim),'ln(z-score power)','fontsize',12,'rotation',90,'horizontalalignment','center')
title('No Stimulus/Speech (z-scored by frequency)','fontweight','normal')
set(gca,'yscale','log','ytick',ft,'yticklabel',ftl);
caxis([minZ,maxZ])

subplot(3,1,2)
pcolor(binz(2:size(mz_zSpeech,3)+1),frx,squeeze(nanmean(mz_zSpeech,2))); shading flat; set(gca,'ydir','normal'); ylabel('Frequency (Hz)'); xlabel('Distance (mm)'); set(gca,'fontsize',14); colorbar;
text(max(xlim)+diff(xlim)/4,mean(ylim),'ln(z-score power)','fontsize',12,'rotation',90,'horizontalalignment','center')
title('Speech (z-scored by frequency) with significant clusters','fontweight','normal')
set(gca,'yscale','log','ytick',ft,'yticklabel',ftl);
% hold on
% for jj = 1:length(clusterSigSpeech)
%    plot(clustXSpeech{jj}(boundarySigSpeech{jj}),clustYSpeech{jj}(boundarySigSpeech{jj}),'k','linewidth',3)   
% end
caxis([minZ,maxZ])

subplot(3,1,3)
pcolor(binz(2:size(mz_zStim,3)+1),frx,squeeze(nanmean(mz_zStim,2))); shading flat; set(gca,'ydir','normal'); ylabel('Frequency (Hz)'); xlabel('Distance (mm)'); set(gca,'fontsize',14); colorbar;
text(max(xlim)+diff(xlim)/4,mean(ylim),'ln(z-score power)','fontsize',12,'rotation',90,'horizontalalignment','center')
title('Stimulus (z-scored by frequency) with significant clusters','fontweight','normal')
set(gca,'yscale','log','ytick',ft,'yticklabel',ftl);
% hold on
% for jj = 1:length(clusterSigStim)
%    plot(clustXStim{jj}(boundarySigStim{jj}),clustYStim{jj}(boundarySigStim{jj}),'k','linewidth',3)   
% end
caxis([minZ,maxZ])
tempFig = gcf;
tempFig.Position = [1000 140 938 1198];

if savePlots
    exportgraphics(tempFig,fullfile(folderFigures,pt,[pt '_zscore_power.png']),'Resolution',600)
    exportgraphics(tempFig,fullfile(folderFigures,pt,[pt '_zscore_power.eps']))
end


%%
avgSpeechBase = squeeze(nanmean(mz_zSpeech,2)) - squeeze(nanmean(mz_zNoST,2));
avgStimBase = squeeze(nanmean(mz_zStim,2)) - squeeze(nanmean(mz_zNoST,2));
avgSpeechStim = squeeze(nanmean(mz_zSpeech,2)) - squeeze(nanmean(mz_zStim,2));
maxZsub = max([avgSpeechBase(:);avgStimBase(:);avgSpeechStim(:)]);
minZsub = min([avgSpeechBase(:);avgStimBase(:);avgSpeechStim(:)]);
%
avgSpeechBaseClust = nan(size(avgSpeechBase));
for jj = 1:length(clusterSigSpeech)
    clusterTemp = clusterSigSpeech{jj};
avgSpeechBaseClust(clusterTemp) = avgSpeechBase(clusterTemp);
end;

avgStimBaseClust = nan(size(avgStimBase));
for jj = 1:length(clusterSigStim)
    clusterTemp = clusterSigStim{jj};
avgStimBaseClust(clusterTemp) = avgStimBase(clusterTemp);
end;

avgSpeechStimClust = nan(size(avgSpeechStim));
for jj = 1:length(clusterSigStimSpeech)
    clusterTemp = clusterSigStimSpeech{jj};
avgSpeechStimClust(clusterTemp) = avgSpeechStim(clusterTemp);
end;
%
figure
subplot(3,2,1)
pcolor(binz(2:size(mz_zNoST,3)+1),frx,avgSpeechBase); shading flat; set(gca,'ydir','normal'); ylabel('Frequency (Hz)'); xlabel('Distance (mm)'); set(gca,'fontsize',14); colorbar;
text(max(xlim)+diff(xlim)/4,mean(ylim),'ln(z-score power)','fontsize',12,'rotation',90,'horizontalalignment','center')
title('Speech - Baseline (z-scored by frequency)','fontweight','normal')
set(gca,'yscale','log','ytick',ft,'yticklabel',ftl);
caxis([minZsub,maxZsub])

subplot(3,2,2)
pcolor(binz(2:size(mz_zNoST,3)+1),frx,avgSpeechBaseClust); shading flat; set(gca,'ydir','normal'); ylabel('Frequency (Hz)'); xlabel('Distance (mm)'); set(gca,'fontsize',14); colorbar;
text(max(xlim)+diff(xlim)/4,mean(ylim),'ln(z-score power)','fontsize',12,'rotation',90,'horizontalalignment','center')
title('Speech - Baseline (z-scored by frequency) Significant differences','fontweight','normal')
set(gca,'yscale','log','ytick',ft,'yticklabel',ftl);
caxis([minZsub,maxZsub])

subplot(3,2,3)
pcolor(binz(2:size(mz_zSpeech,3)+1),frx,avgStimBase); shading flat; set(gca,'ydir','normal'); ylabel('Frequency (Hz)'); xlabel('Distance (mm)'); set(gca,'fontsize',14); colorbar;
text(max(xlim)+diff(xlim)/4,mean(ylim),'ln(z-score power)','fontsize',12,'rotation',90,'horizontalalignment','center')
title('Stimulus - Baseline (z-scored by frequency)','fontweight','normal')
set(gca,'yscale','log','ytick',ft,'yticklabel',ftl);
% hold on
% for jj = 1:length(clusterSigSpeech)
%    plot(clustXSpeech{jj}(boundarySigSpeech{jj}),clustYSpeech{jj}(boundarySigSpeech{jj}),'k','linewidth',3)   
% end
caxis([minZsub,maxZsub])

subplot(3,2,4)
pcolor(binz(2:size(mz_zNoST,3)+1),frx,avgStimBaseClust); shading flat; set(gca,'ydir','normal'); ylabel('Frequency (Hz)'); xlabel('Distance (mm)'); set(gca,'fontsize',14); colorbar;
text(max(xlim)+diff(xlim)/4,mean(ylim),'ln(z-score power)','fontsize',12,'rotation',90,'horizontalalignment','center')
title('Stim - Baseline (z-scored by frequency) Significant differences','fontweight','normal')
set(gca,'yscale','log','ytick',ft,'yticklabel',ftl);
caxis([minZsub,maxZsub])

subplot(3,2,5)
pcolor(binz(2:size(mz_zStim,3)+1),frx,avgSpeechStim); shading flat; set(gca,'ydir','normal'); ylabel('Frequency (Hz)'); xlabel('Distance (mm)'); set(gca,'fontsize',14); colorbar;
text(max(xlim)+diff(xlim)/4,mean(ylim),'ln(z-score power)','fontsize',12,'rotation',90,'horizontalalignment','center')
title('Speech - Stimulus (z-scored by frequency)','fontweight','normal')
set(gca,'yscale','log','ytick',ft,'yticklabel',ftl);
% hold on
% for jj = 1:length(clusterSigStim)
%    plot(clustXStim{jj}(boundarySigStim{jj}),clustYStim{jj}(boundarySigStim{jj}),'k','linewidth',3)   
% end
caxis([minZsub,maxZsub])

subplot(3,2,6)
pcolor(binz(2:size(mz_zNoST,3)+1),frx,avgSpeechStimClust); shading flat; set(gca,'ydir','normal'); ylabel('Frequency (Hz)'); xlabel('Distance (mm)'); set(gca,'fontsize',14); colorbar;
text(max(xlim)+diff(xlim)/4,mean(ylim),'ln(z-score power)','fontsize',12,'rotation',90,'horizontalalignment','center')
title('Speech - Stimulus (z-scored by frequency) Significant Differences','fontweight','normal')
set(gca,'yscale','log','ytick',ft,'yticklabel',ftl);
caxis([minZsub,maxZsub])

tempFig = gcf;
tempFig.Position = [839 109 1408 1229];
if savePlots
    exportgraphics(tempFig,fullfile(folderFigures,pt,[pt '_zscore_power_diff.png']),'Resolution',600)
    exportgraphics(tempFig,fullfile(folderFigures,pt,[pt '_zscore_power_diff.eps']))
end
%%
avgNoSTSTG = nanmean(mz_zNoSTSTG,2);
avgSpeechSTG = nanmean(mz_zSpeechSTG,2);
avgStimSTG = nanmean(mz_zStimSTG,2);
maxZSTG = max([avgNoSTSTG(:);avgSpeechSTG(:);avgStimSTG(:)]);
minZSTG = min([avgNoSTSTG(:);avgSpeechSTG(:);avgStimSTG(:)]);


figure
subplot(3,1,1)
pcolor(binz(2:size(mz_zNoSTSTG,3)+1),frx,squeeze(nanmean(mz_zNoSTSTG,2))); shading flat; set(gca,'ydir','normal'); ylabel('Frequency (Hz)'); xlabel('Distance (mm)'); set(gca,'fontsize',14); colorbar;
text(max(xlim)+diff(xlim)/4,mean(ylim),'ln(z-score power)','fontsize',12,'rotation',90,'horizontalalignment','center')
title('STG - No Stimulus/Speech (z-scored by frequency)','fontweight','normal')
set(gca,'yscale','log','ytick',ft,'yticklabel',ftl);
caxis([minZSTG,maxZSTG])

subplot(3,1,2)
pcolor(binz(2:size(mz_zSpeechSTG,3)+1),frx,squeeze(nanmean(mz_zSpeechSTG,2))); shading flat; set(gca,'ydir','normal'); ylabel('Frequency (Hz)'); xlabel('Distance (mm)'); set(gca,'fontsize',14); colorbar;
text(max(xlim)+diff(xlim)/4,mean(ylim),'ln(z-score power)','fontsize',12,'rotation',90,'horizontalalignment','center')
title('STG - Speech (z-scored by frequency) with significant clusters','fontweight','normal')
set(gca,'yscale','log','ytick',ft,'yticklabel',ftl);
% hold on
% for jj = 1:length(clusterSigSpeechSTG)
%    plot(clustXSpeech{jj}(boundarySigSpeech{jj}),clustYSpeech{jj}(boundarySigSpeech{jj}),'k','linewidth',3)   
% end
caxis([minZSTG,maxZSTG])

subplot(3,1,3)
pcolor(binz(2:size(mz_zStimSTG,3)+1),frx,squeeze(nanmean(mz_zStimSTG,2))); shading flat; set(gca,'ydir','normal'); ylabel('Frequency (Hz)'); xlabel('Distance (mm)'); set(gca,'fontsize',14); colorbar;
text(max(xlim)+diff(xlim)/4,mean(ylim),'ln(z-score power)','fontsize',12,'rotation',90,'horizontalalignment','center')
title('STG  - Stimulus (z-scored by frequency) with significant clusters','fontweight','normal')
set(gca,'yscale','log','ytick',ft,'yticklabel',ftl);
% hold on
% for jj = 1:length(clusterSigStimSTG)
%    plot(clustXStimSTG{jj}(boundarySigStimSTG{jj}),clustYStimSTG{jj}(boundarySigStimSTG{jj}),'k','linewidth',3)   
% end
caxis([minZSTG,maxZSTG])

tempFig = gcf;
tempFig.Position = [1000 140 938 1198];
if savePlots
    exportgraphics(tempFig,fullfile(folderFigures,pt,[pt '_zscore_power_STG.png']),'Resolution',600)
    exportgraphics(tempFig,fullfile(folderFigures,pt,[pt '_zscore_power_STG.eps']))
end
%%
avgSpeechBaseSTG = squeeze(nanmean(mz_zSpeechSTG,2)) - squeeze(nanmean(mz_zNoSTSTG,2));
avgStimBaseSTG = squeeze(nanmean(mz_zStimSTG,2)) - squeeze(nanmean(mz_zNoSTSTG,2));
avgSpeechStimSTG = squeeze(nanmean(mz_zSpeechSTG,2)) - squeeze(nanmean(mz_zStimSTG,2));
maxZsubSTG = max([avgSpeechBaseSTG(:);avgStimBaseSTG(:);avgSpeechStimSTG(:)]);
minZsubSTG = min([avgSpeechBaseSTG(:);avgStimBaseSTG(:);avgSpeechStimSTG(:)]);
%
avgSpeechBaseClustSTG = nan(size(avgSpeechBaseSTG));
for jj = 1:length(clusterSigSpeechSTG)
    clusterTemp = clusterSigSpeechSTG{jj};
avgSpeechBaseClustSTG(clusterTemp) = avgSpeechBaseSTG(clusterTemp);
end;

avgStimBaseClustSTG = nan(size(avgStimBaseSTG));
for jj = 1:length(clusterSigStimSTG)
    clusterTemp = clusterSigStimSTG{jj};
avgStimBaseClustSTG(clusterTemp) = avgStimBaseSTG(clusterTemp);
end;

avgSpeechStimClustSTG = nan(size(avgSpeechStimSTG));
for jj = 1:length(clusterSigStimSpeechSTG)
    clusterTemp = clusterSigStimSpeechSTG{jj};
avgSpeechStimClustSTG(clusterTemp) = avgSpeechStimSTG(clusterTemp);
end;
%
figure
subplot(3,2,1)
pcolor(binz(2:size(mz_zNoSTSTG,3)+1),frx,avgSpeechBaseSTG); shading flat; set(gca,'ydir','normal'); ylabel('Frequency (Hz)'); xlabel('Distance (mm)'); set(gca,'fontsize',14); colorbar;
text(max(xlim)+diff(xlim)/4,mean(ylim),'ln(z-score power)','fontsize',12,'rotation',90,'horizontalalignment','center')
title('STG Speech - Baseline (z-scored by frequency) ','fontweight','normal')
set(gca,'yscale','log','ytick',ft,'yticklabel',ftl);
caxis([minZsubSTG,maxZsubSTG])

subplot(3,2,2)
pcolor(binz(2:size(mz_zNoSTSTG,3)+1),frx,avgSpeechBaseClustSTG); shading flat; set(gca,'ydir','normal'); ylabel('Frequency (Hz)'); xlabel('Distance (mm)'); set(gca,'fontsize',14); colorbar;
text(max(xlim)+diff(xlim)/4,mean(ylim),'ln(z-score power)','fontsize',12,'rotation',90,'horizontalalignment','center')
title('STG Speech - Baseline (z-scored by frequency) Significant differences','fontweight','normal')
set(gca,'yscale','log','ytick',ft,'yticklabel',ftl);
caxis([minZsubSTG,maxZsubSTG])

subplot(3,2,3)
pcolor(binz(2:size(mz_zSpeechSTG,3)+1),frx,avgStimBaseSTG); shading flat; set(gca,'ydir','normal'); ylabel('Frequency (Hz)'); xlabel('Distance (mm)'); set(gca,'fontsize',14); colorbar;
text(max(xlim)+diff(xlim)/4,mean(ylim),'ln(z-score power)','fontsize',12,'rotation',90,'horizontalalignment','center')
title('STG Stimulus - Baseline (z-scored by frequency)','fontweight','normal')
set(gca,'yscale','log','ytick',ft,'yticklabel',ftl);
% hold on
% for jj = 1:length(clusterSigSpeech)
%    plot(clustXSpeech{jj}(boundarySigSpeech{jj}),clustYSpeech{jj}(boundarySigSpeech{jj}),'k','linewidth',3)   
% end
caxis([minZsubSTG,maxZsubSTG])

subplot(3,2,4)
pcolor(binz(2:size(mz_zNoSTSTG,3)+1),frx,avgStimBaseClustSTG); shading flat; set(gca,'ydir','normal'); ylabel('Frequency (Hz)'); xlabel('Distance (mm)'); set(gca,'fontsize',14); colorbar;
text(max(xlim)+diff(xlim)/4,mean(ylim),'ln(z-score power)','fontsize',12,'rotation',90,'horizontalalignment','center')
title('STG Stim - Baseline (z-scored by frequency) Significant differences','fontweight','normal')
set(gca,'yscale','log','ytick',ft,'yticklabel',ftl);
caxis([minZsubSTG,maxZsubSTG])

subplot(3,2,5)
pcolor(binz(2:size(mz_zStimSTG,3)+1),frx,avgSpeechStimSTG); shading flat; set(gca,'ydir','normal'); ylabel('Frequency (Hz)'); xlabel('Distance (mm)'); set(gca,'fontsize',14); colorbar;
text(max(xlim)+diff(xlim)/4,mean(ylim),'ln(z-score power)','fontsize',12,'rotation',90,'horizontalalignment','center')
title('STG Speech - Stimulus (z-scored by frequency)','fontweight','normal')
set(gca,'yscale','log','ytick',ft,'yticklabel',ftl);
% hold on
% for jj = 1:length(clusterSigStim)
%    plot(clustXStim{jj}(boundarySigStim{jj}),clustYStim{jj}(boundarySigStim{jj}),'k','linewidth',3)   
% end
caxis([minZsubSTG,maxZsubSTG])

subplot(3,2,6)
pcolor(binz(2:size(mz_zStimSTG,3)+1),frx,avgSpeechStimClustSTG); shading flat; set(gca,'ydir','normal'); ylabel('Frequency (Hz)'); xlabel('Distance (mm)'); set(gca,'fontsize',14); colorbar;
text(max(xlim)+diff(xlim)/4,mean(ylim),'ln(z-score power)','fontsize',12,'rotation',90,'horizontalalignment','center')
title('STG Speech - Stimulus (z-scored by frequency) Significant Differences','fontweight','normal')
set(gca,'yscale','log','ytick',ft,'yticklabel',ftl);
caxis([minZsubSTG,maxZsubSTG])

tempFig = gcf;
tempFig.Position = [839 109 1408 1229];
if savePlots
    exportgraphics(tempFig,fullfile(folderFigures,pt,[pt '_zscore_power_diff_STG.png']),'Resolution',600)
    exportgraphics(tempFig,fullfile(folderFigures,pt,[pt '_zscore_power_diff_STG.eps']))
end
%%
avgNoSTIFG = nanmean(mz_zNoSTIFG,2);
avgSpeechIFG = nanmean(mz_zSpeechIFG,2);
avgStimIFG = nanmean(mz_zStimIFG,2);
maxZIFG = max([avgNoSTIFG(:);avgSpeechIFG(:);avgStimIFG(:)]);
minZIFG = min([avgNoSTIFG(:);avgSpeechIFG(:);avgStimIFG(:)]);

figure
subplot(3,1,1)
pcolor(binz(2:size(mz_zNoSTIFG,3)+1),frx,squeeze(nanmean(mz_zNoSTIFG,2))); shading flat; set(gca,'ydir','normal'); ylabel('Frequency (Hz)'); xlabel('Distance (mm)'); set(gca,'fontsize',14); colorbar;
text(max(xlim)+diff(xlim)/4,mean(ylim),'ln(z-score power)','fontsize',12,'rotation',90,'horizontalalignment','center')
title('IFG - No Stimulus/Speech (z-scored by frequency)','fontweight','normal')
set(gca,'yscale','log','ytick',ft,'yticklabel',ftl);
caxis([minZIFG,maxZIFG])

subplot(3,1,2)
pcolor(binz(2:size(mz_zSpeechIFG,3)+1),frx,squeeze(nanmean(mz_zSpeechIFG,2))); shading flat; set(gca,'ydir','normal'); ylabel('Frequency (Hz)'); xlabel('Distance (mm)'); set(gca,'fontsize',14); colorbar;
text(max(xlim)+diff(xlim)/4,mean(ylim),'ln(z-score power)','fontsize',12,'rotation',90,'horizontalalignment','center')
title('IFG - Speech (z-scored by frequency) with significant clusters','fontweight','normal')
set(gca,'yscale','log','ytick',ft,'yticklabel',ftl);
% hold on
% for jj = 1:length(clusterSigSpeechIFG)
%    plot(clustXSpeechIFG{jj}(boundarySigSpeechIFG{jj}),clustYSpeechIFG{jj}(boundarySigSpeechIFG{jj}),'k','linewidth',3)   
% end
caxis([minZIFG,maxZIFG])


subplot(3,1,3)
pcolor(binz(2:size(mz_zStimIFG,3)+1),frx,squeeze(nanmean(mz_zStimIFG,2))); shading flat; set(gca,'ydir','normal'); ylabel('Frequency (Hz)'); xlabel('Distance (mm)'); set(gca,'fontsize',14); colorbar;
text(max(xlim)+diff(xlim)/4,mean(ylim),'ln(z-score power)','fontsize',12,'rotation',90,'horizontalalignment','center')
title('IFG - Stimulus (z-scored by frequency) with significant clusters','fontweight','normal')
set(gca,'yscale','log','ytick',ft,'yticklabel',ftl);
% hold on
% for jj = 1:length(clusterSigStimIFG)
%    plot(clustXStimIFG{jj}(boundarySigStimIFG{jj}),clustYStimIFG{jj}(boundarySigStimIFG{jj}),'k','linewidth',3)   
% end
caxis([minZIFG,maxZIFG])

tempFig = gcf;
tempFig.Position = [1000 140 938 1198];
if savePlots
    exportgraphics(tempFig,fullfile(folderFigures,pt,[pt '_zscore_power_IFG.png']),'Resolution',600)
    exportgraphics(tempFig,fullfile(folderFigures,pt,[pt '_zscore_power_IFG.eps']))
end
%%
avgSpeechBaseIFG = squeeze(nanmean(mz_zSpeechIFG,2)) - squeeze(nanmean(mz_zNoSTIFG,2));
avgStimBaseIFG = squeeze(nanmean(mz_zStimIFG,2)) - squeeze(nanmean(mz_zNoSTIFG,2));
avgSpeechStimIFG = squeeze(nanmean(mz_zSpeechIFG,2)) - squeeze(nanmean(mz_zStimIFG,2));
maxZsubIFG = max([avgSpeechBaseIFG(:);avgStimBaseIFG(:);avgSpeechStimIFG(:)]);
minZsubIFG = min([avgSpeechBaseIFG(:);avgStimBaseIFG(:);avgSpeechStimIFG(:)]);
%
avgSpeechBaseClustIFG = nan(size(avgSpeechBaseIFG));
for jj = 1:length(clusterSigSpeechIFG)
    clusterTemp = clusterSigSpeechIFG{jj};
avgSpeechBaseClustIFG(clusterTemp) = avgSpeechBaseIFG(clusterTemp);
end;

avgStimBaseClustIFG = nan(size(avgStimBaseIFG));
for jj = 1:length(clusterSigStimIFG)
    clusterTemp = clusterSigStimIFG{jj};
avgStimBaseClustIFG(clusterTemp) = avgStimBaseIFG(clusterTemp);
end;

avgSpeechStimClustIFG = nan(size(avgSpeechStimIFG));
for jj = 1:length(clusterSigStimSpeechIFG)
    clusterTemp = clusterSigStimSpeechIFG{jj};
avgSpeechStimClustIFG(clusterTemp) = avgSpeechStimIFG(clusterTemp);
end;
%
figure
subplot(3,2,1)
pcolor(binz(2:size(mz_zNoSTIFG,3)+1),frx,avgSpeechBaseIFG); shading flat; set(gca,'ydir','normal'); ylabel('Frequency (Hz)'); xlabel('Distance (mm)'); set(gca,'fontsize',14); colorbar;
text(max(xlim)+diff(xlim)/4,mean(ylim),'ln(z-score power)','fontsize',12,'rotation',90,'horizontalalignment','center')
title('IFG Speech - Baseline (z-scored by frequency) ','fontweight','normal')
set(gca,'yscale','log','ytick',ft,'yticklabel',ftl);
caxis([minZsubIFG,maxZsubIFG])

subplot(3,2,2)
pcolor(binz(2:size(mz_zNoSTIFG,3)+1),frx,avgSpeechBaseClustIFG); shading flat; set(gca,'ydir','normal'); ylabel('Frequency (Hz)'); xlabel('Distance (mm)'); set(gca,'fontsize',14); colorbar;
text(max(xlim)+diff(xlim)/4,mean(ylim),'ln(z-score power)','fontsize',12,'rotation',90,'horizontalalignment','center')
title('IFG Speech - Baseline (z-scored by frequency) Significant differences','fontweight','normal')
set(gca,'yscale','log','ytick',ft,'yticklabel',ftl);
caxis([minZsubIFG,maxZsubIFG])

subplot(3,2,3)
pcolor(binz(2:size(mz_zSpeechIFG,3)+1),frx,avgStimBaseIFG); shading flat; set(gca,'ydir','normal'); ylabel('Frequency (Hz)'); xlabel('Distance (mm)'); set(gca,'fontsize',14); colorbar;
text(max(xlim)+diff(xlim)/4,mean(ylim),'ln(z-score power)','fontsize',12,'rotation',90,'horizontalalignment','center')
title('IFG Stimulus - Baseline (z-scored by frequency)','fontweight','normal')
set(gca,'yscale','log','ytick',ft,'yticklabel',ftl);
% hold on
% for jj = 1:length(clusterSigSpeech)
%    plot(clustXSpeech{jj}(boundarySigSpeech{jj}),clustYSpeech{jj}(boundarySigSpeech{jj}),'k','linewidth',3)   
% end
caxis([minZsubIFG,maxZsubIFG])

subplot(3,2,4)
pcolor(binz(2:size(mz_zNoSTIFG,3)+1),frx,avgStimBaseClustIFG); shading flat; set(gca,'ydir','normal'); ylabel('Frequency (Hz)'); xlabel('Distance (mm)'); set(gca,'fontsize',14); colorbar;
text(max(xlim)+diff(xlim)/4,mean(ylim),'ln(z-score power)','fontsize',12,'rotation',90,'horizontalalignment','center')
title('Stim - Baseline (z-scored by frequency) Significant differences','fontweight','normal')
set(gca,'yscale','log','ytick',ft,'yticklabel',ftl);
caxis([minZsubIFG,maxZsubIFG])

subplot(3,2,5)
pcolor(binz(2:size(mz_zStimIFG,3)+1),frx,avgSpeechStimIFG); shading flat; set(gca,'ydir','normal'); ylabel('Frequency (Hz)'); xlabel('Distance (mm)'); set(gca,'fontsize',14); colorbar;
text(max(xlim)+diff(xlim)/4,mean(ylim),'ln(z-score power)','fontsize',12,'rotation',90,'horizontalalignment','center')
title('Speech - Stimulus (z-scored by frequency)','fontweight','normal')
set(gca,'yscale','log','ytick',ft,'yticklabel',ftl);
% hold on
% for jj = 1:length(clusterSigStim)
%    plot(clustXStim{jj}(boundarySigStim{jj}),clustYStim{jj}(boundarySigStim{jj}),'k','linewidth',3)   
% end
caxis([minZsubIFG,maxZsubIFG])

subplot(3,2,6)
pcolor(binz(2:size(mz_zStimIFG,3)+1),frx,avgSpeechStimClustIFG); shading flat; set(gca,'ydir','normal'); ylabel('Frequency (Hz)'); xlabel('Distance (mm)'); set(gca,'fontsize',14); colorbar;
text(max(xlim)+diff(xlim)/4,mean(ylim),'ln(z-score power)','fontsize',12,'rotation',90,'horizontalalignment','center')
title('Speech - Stimulus (z-scored by frequency) Significant Differences','fontweight','normal')
set(gca,'yscale','log','ytick',ft,'yticklabel',ftl);
caxis([minZsubIFG,maxZsubIFG])
tempFig = gcf;
tempFig.Position = [839 109 1408 1229];
if savePlots
    exportgraphics(tempFig,fullfile(folderFigures,pt,[pt '_zscore_power_diff_IFG.png']),'Resolution',600)
    exportgraphics(tempFig,fullfile(folderFigures,pt,[pt '_zscore_power_diff_IFG.eps']))
end
%%
figure('color','w','position',[277 223 1095 963]); colormap(parula);

% Plot confusion matrix of bipolar distances between all pairs
subplot(8,2,1:2:5); pcolor(Mbpdist); shf; hold on; plot([1 1; 1 256]',[1 256;256 256]','k-','linewidth',1); plot([1 1 1 128 256],[1 128 256 256 256],'k.'); %plot([1 128 256 256 256]/2,[1 1 1 128 256]/2,'k.');
axis equal off; set(gca,'ydir','reverse','xdir','reverse'); colormap(gca,flipud(cmocean('deep'))); colorbar('fontsize',12); title('Bipolar pair distance','fontsize',14,'fontweight','normal')

% Histogram of number of pairs per bin
subplot(8,2,7); histogram(make1d(Mbpdist),0:binsz:85,'facecolor',.5*[1 1 1]); set(gca,'fontsize',12); xlabel('Binned bipolar distance (mm)','fontsize',14);
ylabel('Counts/bin','fontweight','normal'); axis tight; grid on; cb=colorbar; set(cb,'visible','off');

% line plot of same as above, zscored
subplot(2,2,2); hold on; cmo=fliplr(cmocean('thermal',nfrx+round(nfrx)*.1)); colormap(gca,cmo); cb=colorbar; set(cb,'ticks',0:.25:1,'ticklabels',0:50:200,'fontsize',12)
for i=nfrx:-1:1; nns=~isnan(mz(i,:));
    plot(binz(2:size(mz_z,2)+1),mz_z(i,:),'-','color',cmo(i,:),'linewidth',1);
end; grid on; axis tight;
set(gca,'fontsize',14); ylabel('ln(power)'); xlabel('Bipolar distance (mm)','fontsize',14);
text(max(xlim)+diff(xlim)/4,mean(ylim),'Frequency (Hz)','fontsize',12,'rotation',90,'horizontalalignment','center')

subplot(2,2,3);
pcolor(binz(2:size(mz_z,2)+1),frx,mz_z); shading flat; set(gca,'ydir','normal'); ylabel('Frequency (Hz)'); xlabel('Distance (mm)'); set(gca,'fontsize',14); colorbar;
text(max(xlim)+diff(xlim)/4,mean(ylim),'ln(power)','fontsize',12,'rotation',90,'horizontalalignment','center')
title('(z-scored by frequency)','fontweight','normal')
set(gca,'yscale','log','ytick',ft,'yticklabel',ftl);

sp(2,2,4)
[mx,my]=meshgrid(binz(2:size(mz_z,2)+1)+binsz/2,frx);
surf(mx,my,(mz_z)); xlabel('Bipolar distance (mm)','fontsize',14); ylabel('Frequency (Hz)','fontsize',14); zlabel('ln (Power)','fontsize',14)
view(140,25); set(gca,'xdir','reverse','ydir','reverse','ytick',ft,'yticklabel',ftl,'fontsize',14); xlim([binsz/2 85]);  axis tight;
set(gca,'yscale','log');
title('(z-scored by frequency)','fontweight','normal')%set(gca,'yscale','linear'); set(gca,'zscale','log'); set(gca,'zscale','linear')
%ylim([0 30]); %zlim([.035 3]); caxis([.035 3]); %low freqs only
%ylim([30 200]); %zlim([0 0.035]); caxis([0 0.035]); %high freqs only
ylim([0 200]);
