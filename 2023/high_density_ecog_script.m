%%

% total number of electrodes that have spikes, what's the  maximal distance
% between electrode pairs

% same exact analysis except look at 51 - 52 , look at 52 - 53, vs. 51 - 53
% for sub-sampled. so shift it over one - do this later

% show bipolar pairs spectra 4 mm vs 8 mm for baseline - schematic

% look at SVM on includec chnnals that match SVM

% look at 5 , sample without replacement for baselines , do 5 iterations

% which electrodes are the seizure onset ones - Jon will go through it

% for paper - blunt , its a general audience, for localization, epilepsy
% surgery, histogram of electrode numbers per subject vs. total

% each channel, plot distributions for baseline - histogram of windows x
% baseline , color code

%% load data
doPlotIndSubj = true;
doPerm = true;
doInd = true;
plotPermChans = true;
doConn = false;
doSVM = false;
savePlots = true;

% some algorithm = plotting choices, defaults are for raw spikes, for line length
% algorthim the choices are further below
colormapChoice = 'PRgn';
desiredPlotBounds = [-1000 1000];
twoSidedPerm = true;

fs = 512; % sampling rate Hz

data_root = getenv("KLEEN_DATA");
datadir = fullfile(data_root, 'bipolar_expedition');
folderDataBase = '/home/dcaldwell/results'; %%% change to be on the main data drive 


%fileSpikes = '/Volumes/KLEEN_DRIVE/David/Bipolar project/taggedspikes.mat'; %old
fileSpikes = fullfile(datadir, 'taggedspikes_April2022.mat');
folderBaseline = fullfile(datadir,'baseline-high-density-data');
filesFolderBaseline = dir(folderBaseline);
namesFilesBaseline = {filesFolderBaseline(:).name};


load(fileSpikes)
ptsTotal = {'EC129';'EC130';'EC131';'EC133';'EC135';'EC136';'EC137';'EC143';'EC157';'EC162';'EC168';'EC175';'EC181';'EC183';'EC186';'EC187';'EC191';'EC195';'EC196';'EC219';'EC220';'EC221';'EC222'};
% no EC136
ptsTotal = {'EC129';'EC130';'EC131';'EC133';'EC137';'EC143';'EC157';'EC162';'EC168';'EC175';'EC183';'EC186';'EC187';'EC191';'EC195';'EC196';'EC219';'EC220';'EC221';'EC222'};

ptsTotal = {'EC129';'EC130';'EC131';'EC133';'EC137';'EC143';'EC157';'EC162';'EC168';'EC175';'EC183';'EC186';'EC187';'EC191';'EC196';'EC219';'EC220';'EC221';'EC222'};

% exclude 129, 130, 137, 175, 183, 187, 191, 196, 219 because < 50 spikes 
ptsTotal = {'EC131';'EC133';'EC143';'EC157';'EC162';'EC168';'EC186';'EC220';'EC221';'EC222'};


%ptsTotal = {'EC162';'EC168';'EC175';'EC183';'EC186';'EC187';'EC191';'EC196';'EC219';'EC220';'EC221';'EC222'};


%ptsTotal = {'EC196';'EC219';'EC220';'EC221';'EC222'};

%ptsTotal = {'EC219'};

%%
%Spt = patient name
%Straces_allch = spike data
%Selem = 7th column has a 1 if it was in the Intertrial-Interval, 0 if it was during a trial

patientsVetted = ptsTotal;
% look back at 131 - doesnt have any data
%  EC 137 - why are some channels flat
% EC 168 - also why blank

% skipping EC135 for now - no spikes
% skip 136
% EC181 - has no good channels ?
% skip EC195

%folderDataBase = '/Volumes/KLEEN_DRIVE/David/Bipolar project/output';
folderFiguresCell = {fullfile(folderDataBase,'LL20'),fullfile(folderDataBase,'LL40'),fullfile(folderDataBase,'LL100'),fullfile(folderDataBase,'absDer')};
cellAbsDeriv = {false,false,false,true};
cellLL = {true,true,true,false};
cellLLnum = {0.02,0.04,0.1,0};
saveName = {'LL20','LL40','LL100','absDer'};


for index = 1:length(folderFiguresCell)
    
    % visualize spikes to start
    figureMeanSpikes = figure;
    figureMeanReref = figure;
    figureMeanRerefSubSample = figure;
    
    figureBaselineMean = figure;
    figureBaselineMeanReref = figure;
    figureBaselineMeanRerefSubSample = figure;
    
    figureStdSpikes = figure;
    figureStdReref = figure;
    figureStdRerefSubSample = figure;
    
    figureBaselineStd = figure;
    figureBaselineStdReref = figure;
    figureBaselineStdRerefSubSample = figure;
    
    figureTotal = figure;
    tiledlayout(5,5,'TileSpacing','Compact','Padding','Compact');
    hold on;
    figureTotal.Position =  [839 109 1408 1229];
    
    figureTotalSub = figure;
    tiledlayout(5,5,'TileSpacing','Compact','Padding','Compact');
    hold on;
    figureTotalSub.Position =  [839 109 1408 1229];
    
    %permResultsCell = {};
    maxCluster = {};
    permResultsCellConn = {};
    maxClusterConn = {};
    svmCell = {};
    
    % which metric to use for analysis
    useLL = cellLL{index}; % use line length
    useAbsDer = cellAbsDeriv{index};
    folderFigures = folderFiguresCell{index};
    saveNameSpecific = saveName{index};
    
    for jj = 1:length(patientsVetted)
        
        subjName = patientsVetted{jj};
        disp(['Subject ' subjName])
        logicPts = strcmp(subjName,ptsTotal);
        logicVec = find(strcmp(Spt,subjName));
        dataSubj = Straces_allch(logicVec);
        
        % get subject baseline data
        indicesNames = find(contains(namesFilesBaseline,subjName));
        
        totalBaselineDataSubj = [];
        for indicesNums = 1:length(indicesNames)
            baselineData = load(fullfile(folderBaseline,namesFilesBaseline{indicesNames(indicesNums)}));
            spksarti = baselineData.hasspk | baselineData.hasarti;
            excludeInd = spksarti;
            
            % baselineDataV2 = [];
            %baselineCell = baselineData.nonspks_windows(2,:);
            %for numIndices = 1:length(baselineData.nonspks_windows(2,:))
            %baselineDataV2(:,:,numIndices) = baselineCell{numIndices};
            %end
            %baselineDataV2 = permute(baselineDataV2,[2,1,3]);
            
            baselineDataTemp = permute(reshape(cell2mat(baselineData.nonspks_windows(2,~excludeInd)),size(baselineData.nonspks_windows{2,1},1),size(baselineData.nonspks_windows{2,1},2),numel(baselineData.nonspks_windows(2,~excludeInd))),[2,1,3]);
            totalBaselineDataSubj = cat(3,totalBaselineDataSubj,baselineDataTemp);
        end
        % initialize data matrix for each subject
        dataSubjMat = zeros(size(dataSubj{1},1),size(dataSubj{1},2),numel(dataSubj));
        
        for jjj = 1:numel(dataSubj)
            dataSubjMat(:,:,jjj) = dataSubj{jjj};
        end
        
        % make data time x channels x trials
        dataSubjMat = permute(dataSubjMat,[2,1,3]);
        
        % select number of baseline trials to be same as number of trials for
        % each subject
        if size(totalBaselineDataSubj,3) >numel(dataSubj) % account for case where more spikes than baselines
            totalBaselineDataSubj = totalBaselineDataSubj(:,:,randperm(size(totalBaselineDataSubj,3),numel(dataSubj)));
        elseif size(totalBaselineDataSubj,3) <= numel(dataSubj)
            dataSubjMat = dataSubjMat(:,:,randperm(size(dataSubjMat,3),size(totalBaselineDataSubj,3)));
            %    totalBaselineDataSubj = totalBaselineDataSubj; % if same size, just have to match them up
        end
        
        disp(['Number of spike windows ' num2str(numel(dataSubj))])
        disp(['Number of baseline windows ' num2str(size(totalBaselineDataSubj,3))])
        
        %trial indicator
        trialsSubj = Selem(logicVec,7);
        
        t = linspace(-0.5,0.5,size(dataSubjMat,1));
        %dataSubjMat = reshape(cell2mat(dataSubj),size(dataSubj{1},1),numel(dataSubj),size(dataSubj{1},2));
        %dataSubjMat = permute(dataSubjMat,[3,1,2]);
        
        % deal with bad channels
        badChansSubj = badchidx{logicPts};
        goodChansVec = true(size(dataSubj{1},1),1);
        goodChansVec(badChansSubj) = false;
        
        dataSubjMat(:,~goodChansVec,:) = nan;
        
        dataSubjMatGoodOnly = dataSubjMat(:,goodChansVec,:);
        totalBaselineDataSubjGoodOnly = totalBaselineDataSubj(:,goodChansVec,:);
        
        meanSpikesGoodOnly = squeeze(mean(dataSubjMatGoodOnly,3));
        meanSpikes = squeeze(mean(dataSubjMat,3));
        
        stdSpikesGoodOnly = squeeze(std(dataSubjMatGoodOnly,0,3));
        stdSpikes = squeeze(std(dataSubjMat,0,3));
        
        meanBaselineGoodOnly = squeeze(mean(totalBaselineDataSubjGoodOnly,3));
        meanBaseline = squeeze(mean(totalBaselineDataSubj,3));
        
        stdBaselineGoodOnly = squeeze(std(totalBaselineDataSubjGoodOnly,0,3));
        stdBaseline = squeeze(std(totalBaselineDataSubj,0,3));
        
        % rereference
        referencedData = ref2bipolar_reg(dataSubjMat,subjName);
        referencedDataSubSample = ref2bipolar_subsample(dataSubjMat,subjName);
        
        referencedDataGoodOnly = referencedData(:,goodChansVec,:);
        referencedDataSubSampleGoodOnly = referencedDataSubSample(:,goodChansVec,:);
        
        meanSpikesReferencedGoodOnly = squeeze(mean(referencedDataGoodOnly,3));
        meanSpikesReferencedSubSampleGoodOnly = squeeze(mean(referencedDataSubSampleGoodOnly,3));
        
        stdSpikesReferencedGoodOnly = squeeze(std(referencedDataGoodOnly,0,3));
        stdSpikesReferencedSubSampleGoodOnly = squeeze(std(referencedDataSubSampleGoodOnly,0,3));
        
        % rereference baseline
        referencedDataBaseline = ref2bipolar_reg(totalBaselineDataSubj,subjName);
        referencedDataBaselineSubSample = ref2bipolar_subsample(totalBaselineDataSubj,subjName);
        
        referencedDataBaselineGoodOnly = referencedDataBaseline(:,goodChansVec,:);
        referencedDataBaselineSubSampleGoodOnly = referencedDataBaselineSubSample(:,goodChansVec,:);
        
        meanSpikesReferencedBaselineGoodOnly = squeeze(mean(referencedDataBaselineGoodOnly,3));
        meanSpikesReferencedBaselineSubSampleGoodOnly = squeeze(mean(referencedDataBaselineSubSampleGoodOnly,3));
        
        stdSpikesReferencedBaselineGoodOnly = squeeze(std(referencedDataBaselineGoodOnly,0,3));
        stdSpikesReferencedBaselineSubSampleGoodOnly = squeeze(std(referencedDataBaselineSubSampleGoodOnly,0,3));
        
        %
        notNanCombined = ~isnan(referencedDataSubSampleGoodOnly) & ~isnan(referencedDataGoodOnly);
        notNanIndsCombined = find(squeeze(notNanCombined(256,:,1)));
        
        %
        if useLL
            llw = cellLLnum{index};
            % line length processing
            referencedDataLL = line_length_alg(referencedData,fs,llw);
            referencedDataSubSampleLL = line_length_alg(referencedDataSubSample,fs,llw);
            referencedDataGoodOnlyLL = line_length_alg(referencedDataGoodOnly,fs,llw);
            referencedDataSubSampleGoodOnlyLL = line_length_alg(referencedDataSubSampleGoodOnly,fs,llw);
            
            referencedDataBaselineLL= line_length_alg(referencedDataBaseline,fs,llw);
            referencedDataBaselineSubSampleLL= line_length_alg(referencedDataBaselineSubSample,fs,llw);
            referencedDataBaselineGoodOnlyLL= line_length_alg(referencedDataBaselineGoodOnly,fs,llw);
            referencedDataBaselineSubSampleGoodOnlyLL= line_length_alg(referencedDataBaselineSubSampleGoodOnly,fs,llw);
            
            referencedData = referencedDataLL;
            referencedDataGoodOnly = referencedDataGoodOnlyLL;
            referencedDataSubSample = referencedDataSubSampleLL;
            referencedDataSubSampleGoodOnly = referencedDataSubSampleGoodOnlyLL;
            
            referencedDataBaseline = referencedDataBaselineLL;
            referencedDataBaselineSubSample = referencedDataBaselineSubSampleLL;
            referencedDataBaselineGoodOnly = referencedDataBaselineGoodOnlyLL;
            referencedDataBaselineSubSampleGoodOnly = referencedDataBaselineSubSampleGoodOnlyLL;
            
            meanSpikesReferencedGoodOnly = squeeze(mean(referencedDataGoodOnly,3));
            meanSpikesReferencedSubSampleGoodOnly = squeeze(mean(referencedDataSubSampleGoodOnly,3));
            
            % algorithm and plotting choices for line length
            colormapChoice = 'PuRd';
            desiredPlotBounds = [0 200];
            twoSidedPerm = false;
            
        end
        
        if useAbsDer
            % absolute value derivative
            referencedDataAbsDer = abs_deriv(referencedData);
            referencedDataSubSampleAbsDer = abs_deriv(referencedDataSubSample);
            referencedDataGoodOnlyAbsDer = abs_deriv(referencedDataGoodOnly);
            referencedDataSubSampleGoodOnlyAbsDer = abs_deriv(referencedDataSubSampleGoodOnly);
            
            referencedDataBaselineAbsDer = abs_deriv(referencedDataBaseline);
            referencedDataBaselineSubSampleAbsDer = abs_deriv(referencedDataBaselineSubSample);
            referencedDataBaselineGoodOnlyAbsDer = abs_deriv(referencedDataBaselineGoodOnly);
            referencedDataBaselineSubSampleGoodOnlyAbsDer = abs_deriv(referencedDataBaselineSubSampleGoodOnly);
            
            referencedData = referencedDataAbsDer;
            referencedDataGoodOnly = referencedDataGoodOnlyAbsDer;
            referencedDataSubSample = referencedDataSubSampleAbsDer;
            referencedDataSubSampleGoodOnly = referencedDataSubSampleGoodOnlyAbsDer;
            
            
            referencedDataBaseline = referencedDataBaselineAbsDer;
            referencedDataBaselineSubSample = referencedDataBaselineSubSampleAbsDer;
            referencedDataBaselineGoodOnly = referencedDataBaselineGoodOnlyAbsDer;
            referencedDataBaselineSubSampleGoodOnly = referencedDataBaselineSubSampleGoodOnlyAbsDer;
            
            meanSpikesReferencedGoodOnly = squeeze(mean(referencedDataGoodOnly,3));
            meanSpikesReferencedSubSampleGoodOnly = squeeze(mean(referencedDataSubSampleGoodOnly,3));
            
            % algorithm and plotting choices for line length
            colormapChoice = 'PuRd';
            desiredPlotBounds = [0 40];
            twoSidedPerm = false;
            
        end
        
        disp('Data Processed')
        
        %% plotting
        % mean rereferenced data
        figure(figureMeanSpikes)
        subplot(5,5,jj)
        imagesc(meanSpikesGoodOnly(1:end-1,notNanIndsCombined)')
        xt = get(gca, 'XTick');                                             % Original 'XTick' Values
        xtlbl = linspace(-0.5,0.5, numel(xt));                     % New 'XTickLabel' Vector
        set(gca, 'XTick',xt, 'XTickLabel',xtlbl)   % Label Ticks
        title(['Subject ' subjName])
        % caxis([-max(abs(meanSpikesGoodOnly(:))) max(abs(meanSpikesGoodOnly(:)))])
        caxis(desiredPlotBounds)
        colormap(brewermap([],colormapChoice))
        
        figure(figureMeanRerefSubSample)
        subplot(5,5,jj)
        imagesc(meanSpikesReferencedSubSampleGoodOnly(1:end-1,notNanIndsCombined)')
        xt = get(gca, 'XTick');                                             % Original 'XTick' Values
        xtlbl = linspace(-0.5,0.5, numel(xt));                     % New 'XTickLabel' Vector
        set(gca, 'XTick',xt, 'XTickLabel',xtlbl)   % Label Ticks
        title(['Subject ' subjName])
        %  caxis([-max(abs(meanSpikesReferencedSubSampleGoodOnly(:))) max(abs(meanSpikesReferencedSubSampleGoodOnly(:)))])
        caxis(desiredPlotBounds)
        colormap(brewermap([],colormapChoice))
        
        figure(figureMeanReref)
        subplot(5,5,jj)
        imagesc(meanSpikesReferencedGoodOnly(1:end-1,notNanIndsCombined)')
        xt = get(gca, 'XTick');                                             % Original 'XTick' Values
        xtlbl = linspace(-0.5,0.5, numel(xt));                     % New 'XTickLabel' Vector
        set(gca, 'XTick',xt, 'XTickLabel',xtlbl)   % Label Ticks
        title(['Subject ' subjName])
        %caxis([-max(abs(meanSpikesReferencedGoodOnly(:))) max(abs(meanSpikesReferencedGoodOnly(:)))])
        caxis(desiredPlotBounds)
        colormap(brewermap([],colormapChoice))
        
        % mean baseline data
        figure(figureBaselineMean)
        subplot(5,5,jj)
        imagesc(meanBaselineGoodOnly(:,notNanIndsCombined)')
        xt = get(gca, 'XTick');                                             % Original 'XTick' Values
        xtlbl = linspace(-0.5,0.5, numel(xt));                     % New 'XTickLabel' Vector
        set(gca, 'XTick',xt, 'XTickLabel',xtlbl)   % Label Ticks
        title(['Subject ' subjName])
        %caxis([-max(abs(meanBaselineGoodOnly(:))) max(abs(meanBaselineGoodOnly(:)))])
        caxis(desiredPlotBounds)
        colormap(brewermap([],colormapChoice))
        
        figure(figureBaselineMeanRerefSubSample)
        subplot(5,5,jj)
        imagesc(meanSpikesReferencedBaselineSubSampleGoodOnly(:,notNanIndsCombined)')
        xt = get(gca, 'XTick');                                             % Original 'XTick' Values
        xtlbl = linspace(-0.5,0.5, numel(xt));                     % New 'XTickLabel' Vector
        set(gca, 'XTick',xt, 'XTickLabel',xtlbl)   % Label Ticks
        title(['Subject ' subjName])
        %caxis([-max(abs(meanSpikesReferencedBaselineSubSampleGoodOnly(:))) max(abs(meanSpikesReferencedBaselineSubSampleGoodOnly(:)))])
        caxis(desiredPlotBounds)
        colormap(brewermap([],colormapChoice))
        
        figure(figureBaselineMeanReref)
        subplot(5,5,jj)
        imagesc(meanSpikesReferencedBaselineGoodOnly(:,notNanIndsCombined)')
        xt = get(gca, 'XTick');                                             % Original 'XTick' Values
        xtlbl = linspace(-0.5,0.5, numel(xt));                     % New 'XTickLabel' Vector
        set(gca, 'XTick',xt, 'XTickLabel',xtlbl)   % Label Ticks
        title(['Subject ' subjName])
        %caxis([-max(abs(meanSpikesReferencedBaselineGoodOnly(:))) max(abs(meanSpikesReferencedBaselineGoodOnly(:)))])
        caxis(desiredPlotBounds)
        colormap(brewermap([],colormapChoice))
        
        % std rereferenced data
        figure(figureStdSpikes)
        subplot(5,5,jj)
        imagesc(stdSpikesGoodOnly(1:end-1,notNanIndsCombined)')
        xt = get(gca, 'XTick');                                             % Original 'XTick' Values
        xtlbl = linspace(-0.5,0.5, numel(xt));                     % New 'XTickLabel' Vector
        set(gca, 'XTick',xt, 'XTickLabel',xtlbl)   % Label Ticks
        title(['Subject ' subjName])
        % caxis([-max(abs(stdSpikesGoodOnly(:))) max(abs(stdSpikesGoodOnly(:)))])
        caxis([-2000 2000])
        colormap(brewermap([],colormapChoice))
        
        figure(figureStdRerefSubSample)
        subplot(5,5,jj)
        imagesc(stdSpikesReferencedSubSampleGoodOnly(1:end-1,notNanIndsCombined)')
        xt = get(gca, 'XTick');                                             % Original 'XTick' Values
        xtlbl = linspace(-0.5,0.5, numel(xt));                     % New 'XTickLabel' Vector
        set(gca, 'XTick',xt, 'XTickLabel',xtlbl)   % Label Ticks
        
        title(['Subject ' subjName])
        %caxis([-max(abs(stdSpikesReferencedSubSampleGoodOnly(:))) max(abs(stdSpikesReferencedSubSampleGoodOnly(:)))])
        caxis([-2000 2000])
        colormap(brewermap([],colormapChoice))
        
        figure(figureStdReref)
        subplot(5,5,jj)
        imagesc(stdSpikesReferencedGoodOnly(1:end-1,notNanIndsCombined)')
        xt = get(gca, 'XTick');                                             % Original 'XTick' Values
        xtlbl = linspace(-0.5,0.5, numel(xt));                     % New 'XTickLabel' Vector
        set(gca, 'XTick',xt, 'XTickLabel',xtlbl)   % Label Ticks
        title(['Subject ' subjName])
        % caxis([-max(abs(stdSpikesReferencedGoodOnly(:))) max(abs(stdSpikesReferencedGoodOnly(:)))])
        caxis([-2000 2000])
        colormap(brewermap([],colormapChoice))
        
        % std baseline data
        figure(figureBaselineStd)
        subplot(5,5,jj)
        imagesc(stdBaselineGoodOnly(:,notNanIndsCombined)')
        xt = get(gca, 'XTick');                                             % Original 'XTick' Values
        xtlbl = linspace(-0.5,0.5, numel(xt));                     % New 'XTickLabel' Vector
        set(gca, 'XTick',xt, 'XTickLabel',xtlbl)   % Label Ticks
        title(['Subject ' subjName])
        caxis([-max(abs(stdBaselineGoodOnly(:))) max(abs(stdBaselineGoodOnly(:)))])
        caxis([-2000 2000])
        colormap(brewermap([],colormapChoice))
        
        figure(figureBaselineStdRerefSubSample)
        subplot(5,5,jj)
        imagesc(stdSpikesReferencedBaselineSubSampleGoodOnly(:,notNanIndsCombined)')
        xt = get(gca, 'XTick');                                             % Original 'XTick' Values
        xtlbl = linspace(-0.5,0.5, numel(xt));                     % New 'XTickLabel' Vector
        set(gca, 'XTick',xt, 'XTickLabel',xtlbl)   % Label Ticks
        title(['Subject ' subjName])
        caxis([-max(abs(stdSpikesReferencedBaselineSubSampleGoodOnly(:))) max(abs(stdSpikesReferencedBaselineSubSampleGoodOnly(:)))])
        caxis([-2000 2000])
        colormap(brewermap([],colormapChoice))
        
        figure(figureBaselineStdReref)
        subplot(5,5,jj)
        imagesc(stdSpikesReferencedBaselineGoodOnly(:,notNanIndsCombined)')
        xt = get(gca, 'XTick');                                             % Original 'XTick' Values
        xtlbl = linspace(-0.5,0.5, numel(xt));                     % New 'XTickLabel' Vector
        set(gca, 'XTick',xt, 'XTickLabel',xtlbl)   % Label Ticks
        title(['Subject ' subjName])
        %  caxis([-max(abs(stdSpikesReferencedBaselineGoodOnly(:))) max(abs(stdSpikesReferencedBaselineGoodOnly(:)))])
        caxis([-2000 2000])
        colormap(brewermap([],colormapChoice))
        
        if doPlotIndSubj
            figure
            imagesc(meanSpikesReferencedSubSampleGoodOnly(:,notNanIndsCombined)')
            xt = get(gca, 'XTick');                                             % Original 'XTick' Values
            xtlbl = linspace(-0.5,0.5, numel(xt));                     % New 'XTickLabel' Vector
            set(gca, 'XTick',xt, 'XTickLabel',xtlbl)   % Label Ticks
            xlabel('Time (s)')
            ylabel('Channel')
            title([patientsVetted{jj} ' Average spikes - 8 mm rereference'])
            colormap(brewermap([],colormapChoice))
            
            figure
            imagesc(meanSpikesReferencedGoodOnly(:,notNanIndsCombined)')
            xt = get(gca, 'XTick');                                             % Original 'XTick' Values
            xtlbl = linspace(-0.5,0.5, numel(xt));                     % New 'XTickLabel' Vector
            set(gca, 'XTick',xt, 'XTickLabel',xtlbl)   % Label Ticks
            xlabel('Time (s)')
            ylabel('Channel')
            title([patientsVetted{jj} ' Average spikes - 4 mm rereference'])
            colormap(brewermap([],colormapChoice))
            
        end
        %% permutation testing
        if doPerm
            disp(['Doing individual channelwise permutation test'])
            notNan = ~isnan(referencedDataGoodOnly);
            notNanInds = find(squeeze(notNan(256,:,1)));
            notNanSub = ~isnan(referencedDataSubSampleGoodOnly);
            notNanSubInds = find(squeeze(notNanSub(256,:,1)));
            % permutation test on a by channel basis
            
            if doInd
                maxCluster{jj}.subj = subjName;
                tempPvalues = nan(size(referencedDataGoodOnly,2),1);
                % permutation test for background cluster based
                permResultsCell{jj}.subj = subjName;
                
                numGoodChansFull = length(notNanInds);
                tempCellclust = cell(numGoodChansFull,1);
                tempCellpVal = cell(numGoodChansFull,1);
                tempCellpValAdjust = cell(numGoodChansFull,1);
                tempPvaluesInd = nan(numGoodChansFull,1);
                
                numGoodChansFullSub = length(notNanSubInds);
                tempCellclustSub = cell(numGoodChansFullSub,1);
                tempCellpValSub = cell(numGoodChansFullSub,1);
                tempCellpValAdjustSub = cell(numGoodChansFullSub,1);
                tempPvaluesIndSub = nan(numGoodChansFullSub,1);
                
                for chanGoodInd = 1:numGoodChansFull
                    goodChan = notNanInds(chanGoodInd);
                    [clusters, pValues, ~, ~] = permutest(squeeze(referencedDataGoodOnly(1:end-1,goodChan,:)),squeeze(referencedDataBaselineGoodOnly(:,goodChan,:)),false,[],[],twoSidedPerm);
                    
                    tempCellclust{chanGoodInd} = clusters;
                    tempCellpVal{chanGoodInd} = pValues;
                    
                    % build p value matrix
                    tempPvaluesInd(chanGoodInd) = pValues(1);
                    
                    % channelwise FDR
                    [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(pValues,0.05,'dep');
                    tempCellpValAdjust{chanGoodInd} = adj_p;
                    
                    %                                 % extract clusters and tsums
                    %                                 if ~isempty(permResultsCell{jj}.clusters{goodChan}{1})
                    %                                     maxCluster{jj}.clusterMax{goodChan} = clusters{1};
                    %                                     maxCluster{jj}.tSumsMax{goodChan} = tSums(1);
                    %                                     [row,col,zdir] = ind2sub(size(referencedDataGoodOnly(1:end-1,:,:)),clusters{1});
                    %                                     maxCluster{jj}.clusterMaxMatrixInds{goodChan} = [row col zdir];
                    %                                     [boundInd] = boundary(row,col) ;
                    %                                     maxCluster{jj}.boundInd{goodChan} = boundInd;
                    %                                 end
                    disp(['regular sampled done for channel ' num2str(goodChan)])
                end
                
                for chanGoodInd = 1:numGoodChansFull
                    goodChan = notNanInds(chanGoodInd);
                    permResultsCell{jj}.clusters{goodChan} = tempCellclust{chanGoodInd};
                    permResultsCell{jj}.pValues{goodChan} = tempCellpVal{chanGoodInd};
                    permResultsCell{jj}.adjpValues{goodChan} = tempCellpValAdjust{chanGoodInd};
                    tempPvalues(goodChan) = tempPvaluesInd(chanGoodInd);
                end
                
                % do FDR correction across all good channels number 1 p-value
                [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(tempPvalues(notNanIndsCombined),0.05,'dep');
                listChans = [];
                meanWidthList = [];
                for chanGoodInd = 1:length(notNanIndsCombined)
                    goodChan = notNanIndsCombined(chanGoodInd);
                    permResultsCell{jj}.adjpValuesTotal{goodChan} = adj_p(chanGoodInd);
                    permResultsCell{jj}.clusterWidth{goodChan} = max(permResultsCell{jj}.clusters{goodChan}{1})- min(permResultsCell{jj}.clusters{goodChan}{1});
                    if adj_p(chanGoodInd) <= 0.05
                        listChans= [listChans goodChan];
                        meanWidthList = [meanWidthList permResultsCell{jj}.clusterWidth{goodChan}];
                    end
                end
                
                permResultsCell{jj}.numSigChannels = length(listChans);
                permResultsCell{jj}.sigChannels = listChans;
                permResultsCell{jj}.meanWidthList = meanWidthList;
                permResultsCell{jj}.meanWidth = mean(meanWidthList);
                
                tempPvaluesSub = nan(size(referencedDataGoodOnly,2),1);
                
                for chanGoodIndSub = 1:numGoodChansFullSub
                    goodChanSub = notNanSubInds(chanGoodIndSub);
                    [clustersSub, pValuesSub,~, ~]  = permutest(squeeze(referencedDataSubSampleGoodOnly(1:end-1,goodChanSub,:)),squeeze(referencedDataBaselineSubSampleGoodOnly(:,goodChanSub,:)),false,[],[],twoSidedPerm);
                    
                    tempCellclustSub{chanGoodIndSub} = clustersSub;
                    tempCellpValSub{chanGoodIndSub} = pValuesSub;
                    
                    % build p value matrix
                    tempPvaluesIndSub(chanGoodIndSub) = pValuesSub(1);
                    
                    % channelwise FDR
                    [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(pValuesSub,0.05,'dep');
                    tempCellpValAdjustSub{chanGoodIndSub} = adj_p;
                    
                    disp(['subsampled done for channel ' num2str(goodChanSub)])
                    
                    %                                 if ~isempty(permResultsCell{jj}.clustersSub{goodChanSub}{1})
                    %
                    %                                     maxCluster{jj}.clusterMaxSub{goodChan} = clustersSub{1};
                    %                                     maxCluster{jj}.tSumsMaxSub{goodChan} = tSumsSubConn(1);
                    %                                     [row,col,zdir] = ind2sub(size(referencedDataSubSampleGoodOnly(1:end-1,:,:)),clustersSub{1});
                    %                                     maxClusterConn{jj}.clusterMaxMatrixIndsSub{goodChan} = [row col zdir];
                    %                                     [boundInd] = boundary(row,col) ;
                    %                                     maxClusterConn{jj}.boundIndSub{goodChan} = boundInd;
                    %                                 end
                    
                end
                
                for chanGoodIndSub = 1:numGoodChansFullSub
                    goodChanSub = notNanSubInds(chanGoodIndSub);
                    permResultsCell{jj}.clustersSub{goodChanSub} = tempCellclustSub{chanGoodIndSub};
                    permResultsCell{jj}.pValuesSub{goodChanSub} = tempCellpValSub{chanGoodIndSub};
                    permResultsCell{jj}.adjpValuesSub{goodChanSub} = tempCellpValAdjustSub{chanGoodIndSub};
                    tempPvaluesSub(goodChanSub) = tempPvaluesIndSub(chanGoodIndSub);
                end
                
                % do FDR correction across all good channels number 1 p-value
                [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(tempPvaluesSub(notNanIndsCombined),0.05,'dep');
                listChansSub = [];
                meanWidthSubList = [];
                for chanGoodIndSub = 1:length(notNanIndsCombined)
                    goodChanSub = notNanIndsCombined(chanGoodIndSub);
                    permResultsCell{jj}.adjpValuesTotalSub{goodChanSub} = adj_p(chanGoodIndSub);
                    permResultsCell{jj}.clusterWidthSub{goodChanSub} = max(permResultsCell{jj}.clustersSub{goodChanSub}{1})- min(permResultsCell{jj}.clustersSub{goodChanSub}{1});
                    if adj_p(chanGoodIndSub) <= 0.05
                        listChansSub = [listChansSub goodChanSub];
                        meanWidthSubList = [meanWidthSubList permResultsCell{jj}.clusterWidthSub{goodChanSub}];
                    end
                end
                
                permResultsCell{jj}.numSigChannelsSub = length(listChansSub);
                permResultsCell{jj}.sigChannelsSub = listChansSub;
                permResultsCell{jj}.meanWidthListSub = meanWidthSubList;
                permResultsCell{jj}.meanWidthSub = mean(meanWidthSubList);
                
                disp(['Individual Channelwise Permutation test done for Subject ' subjName])
                %%
                dataIntTemp = meanSpikesReferencedGoodOnly(1:end-1,notNanIndsCombined);
                
                figure(figureTotal)
                nexttile
                title(patientsVetted{jj})
                
                figure(figureTotalSub)
                nexttile
                title(patientsVetted{jj})
                
                if plotPermChans & ~isempty(dataIntTemp) & sum(~isnan(dataIntTemp(:)))
                   
                    desiredPlotBoundsSubj = [0 max(dataIntTemp(:))];
                    
                    numPlotVec = length(notNanIndsCombined);
                    tPlot = t(1:end-1);
                    
                    figure
                    imagesc([min(tPlot) max(tPlot)],[1,numPlotVec],meanSpikesReferencedGoodOnly(1:end-1,notNanIndsCombined)')
                    %                xt = get(gca, 'XTick');                                             % Original 'XTick' Values
                    %                xtnew = linspace(min(xt), max(xt), 5);                              % New 'XTick' Values
                    %                xtlbl = linspace(-0.5, 0.5, numel(xtnew));                  % New 'XTickLabel' Vector
                    %                set(gca, 'XTick',xtnew, 'XTickLabel',xtlbl)                         % Label Ticks
                    xlabel('Time (s)')
                    ylabel('Channel')
                    colormap(brewermap([],colormapChoice))
                    caxis(desiredPlotBoundsSubj)
                    hold on
                    for chanGoodInd = 1:length(notNanIndsCombined)
                        goodChan = notNanIndsCombined(chanGoodInd);
                        if ~isempty(permResultsCell{jj}.clusters{goodChan})
                            if ((~isempty(permResultsCell{jj}.clusters{goodChan}{1})) & ( permResultsCell{jj}.adjpValuesTotal{goodChan}(1)<0.05))
                                rectangle('position', [tPlot(min(permResultsCell{jj}.clusters{goodChan}{1})), chanGoodInd-0.5,tPlot(max(permResultsCell{jj}.clusters{goodChan}{1}))- tPlot(min(permResultsCell{jj}.clusters{goodChan}{1})), 1], 'edgecolor', 'k','linewidth',2)
                            end
                        end
                    end
                    currentFig = gcf;
                    currentFig.Position = [839 109 1408 1229];
                    title(['Subject ' subjName ' Average spikes - 4mm rereference - FDR'])
                    set(gca,'FontSize',20)
                    if savePlots
                        exportgraphics(currentFig,fullfile(folderFigures,['fdr_' subjName '_perm.png']),'Resolution',600)
                        exportgraphics(currentFig,fullfile(folderFigures,['fdr_' subjName '_perm.eps']))
                    end
                    
                    
                    figure(figureTotal)
                    hold on
                    imagesc([min(tPlot) max(tPlot)],[1,numPlotVec],meanSpikesReferencedGoodOnly(1:end-1,notNanIndsCombined)')
                    colormap(brewermap([],colormapChoice))
                    caxis(desiredPlotBoundsSubj)
                    hold on
                    for chanGoodInd = 1:length(notNanIndsCombined)
                        goodChan = notNanIndsCombined(chanGoodInd);
                        if ~isempty(permResultsCell{jj}.clusters{goodChan})
                            if ((~isempty(permResultsCell{jj}.clusters{goodChan}{1})) & ( permResultsCell{jj}.adjpValuesTotal{goodChan}(1)<0.05))
                                rectangle('position', [tPlot(min(permResultsCell{jj}.clusters{goodChan}{1})), chanGoodInd-0.5,tPlot(max(permResultsCell{jj}.clusters{goodChan}{1}))- tPlot(min(permResultsCell{jj}.clusters{goodChan}{1})), 1], 'edgecolor', 'k','linewidth',2)
                            end
                        end
                    end
                    title(patientsVetted{jj});
                    
                    
                    figure
                    imagesc([-0.5 0.5],[1,numPlotVec],meanSpikesReferencedGoodOnly(1:end-1,notNanIndsCombined)')
                    %   xticks([-0.5 -0.25 0 0.25 0.5])
                    %   xt = get(gca, 'XTick');                                             % Original 'XTick' Values
                    %  xtlbl = linspace(-0.5,0.5, numel(xt)-1);                     % New 'XTickLabel' Vector
                    % set(gca, 'XTick',xt, 'XTickLabel',xtlbl)   % Label Ticks
                    xlabel('Time (s)')
                    ylabel('Channel')
                    colormap(brewermap([],colormapChoice))
                    caxis(desiredPlotBoundsSubj)
                    hold on
                    for chanGoodInd = 1:length(notNanIndsCombined)
                        goodChan = notNanIndsCombined(chanGoodInd);
                        if ~isempty(permResultsCell{jj}.clusters{goodChan})
                            if ((~isempty(permResultsCell{jj}.clusters{goodChan}{1})) & (permResultsCell{jj}.pValues{goodChan}(1) <0.05))
                                rectangle('position', [tPlot(min(permResultsCell{jj}.clusters{goodChan}{1})), chanGoodInd-0.5,tPlot(max(permResultsCell{jj}.clusters{goodChan}{1}))-tPlot(min(permResultsCell{jj}.clusters{goodChan}{1})), 1], 'edgecolor', 'k','linewidth',2)
                            end
                        end
                    end
                    currentFig = gcf;
                    currentFig.Position = [839 109 1408 1229];
                    title(['Subject ' subjName  ' Average spikes - 4 mm rereference'])
                    set(gca,'FontSize',20)
                    
                    if savePlots
                        exportgraphics(currentFig,fullfile(folderFigures,[subjName '_perm.png']),'Resolution',600)
                        exportgraphics(currentFig,fullfile(folderFigures,[subjName '_perm.eps']))
                    end
                    
                    
                    figure
                    imagesc([-0.5 0.5],[1,numPlotVec],meanSpikesReferencedSubSampleGoodOnly(1:end-1,notNanIndsCombined)')
                    %  xticks([-0.5 -0.25 0 0.25 0.5])
                    %   xt = get(gca, 'XTick');                                             % Original 'XTick' Values
                    %  xtlbl = linspace(-0.5,0.5, numel(xt)-1);                     % New 'XTickLabel' Vector
                    % set(gca, 'XTick',xt, 'XTickLabel',xtlbl)   % Label Ticks
                    xlabel('Time (s)')
                    ylabel('Channel')
                    colormap(brewermap([],colormapChoice))
                    caxis(desiredPlotBoundsSubj)
                    hold on
                    for chanGoodIndSub = 1:length(notNanIndsCombined)
                        goodChanSub = notNanIndsCombined(chanGoodIndSub);
                        if ~isempty(permResultsCell{jj}.clustersSub{goodChanSub})
                            if ((~isempty(permResultsCell{jj}.clustersSub{goodChanSub}{1})) & (permResultsCell{jj}.adjpValuesTotalSub{goodChanSub}(1)<0.05))
                                rectangle('position', [tPlot(min(permResultsCell{jj}.clustersSub{goodChanSub}{1})), chanGoodIndSub-0.5,tPlot(max(permResultsCell{jj}.clustersSub{goodChanSub}{1}))-tPlot(min(permResultsCell{jj}.clustersSub{goodChanSub}{1})), 1], 'edgecolor', 'k','linewidth',2)
                            end
                        end
                    end
                    
                    currentFig = gcf;
                    currentFig.Position = [839 109 1408 1229];
                    
                    title(['Subject ' subjName ' Average spikes - 8 mm rereference - FDR'])
                    set(gca,'FontSize',20)
                    
                    if savePlots
                        exportgraphics(currentFig,fullfile(folderFigures,['fdr_subsample' subjName '_perm.png']),'Resolution',600)
                        exportgraphics(currentFig,fullfile(folderFigures,['fdr_subsample' subjName '_perm.eps']))
                    end
                    
                    
                    figure(figureTotalSub)
                    imagesc([-0.5 0.5],[1,numPlotVec],meanSpikesReferencedSubSampleGoodOnly(1:end-1,notNanIndsCombined)')
                    colormap(brewermap([],colormapChoice))
                    caxis(desiredPlotBoundsSubj)
                    hold on
                    for chanGoodIndSub = 1:length(notNanIndsCombined)
                        goodChanSub = notNanIndsCombined(chanGoodIndSub);
                        if ~isempty(permResultsCell{jj}.clustersSub{goodChanSub})
                            if ((~isempty(permResultsCell{jj}.clustersSub{goodChanSub}{1})) & (permResultsCell{jj}.adjpValuesTotalSub{goodChanSub}(1)<0.05))
                                rectangle('position', [tPlot(min(permResultsCell{jj}.clustersSub{goodChanSub}{1})), chanGoodIndSub-0.5,tPlot(max(permResultsCell{jj}.clustersSub{goodChanSub}{1}))-tPlot(min(permResultsCell{jj}.clustersSub{goodChanSub}{1})), 1], 'edgecolor', 'k','linewidth',2)
                            end
                        end
                    end
                    title(patientsVetted{jj});
                    
                    
                    
                    figure
                    imagesc([-0.5 0.5],[1,numPlotVec],meanSpikesReferencedSubSampleGoodOnly(1:end-1,notNanIndsCombined)')
                    %  xticks([-0.5 -0.25 0 0.25 0.5])
                    %   xt = get(gca, 'XTick');                                             % Original 'XTick' Values
                    %  xtlbl = linspace(-0.5,0.5, numel(xt)-1);                     % New 'XTickLabel' Vector
                    % set(gca, 'XTick',xt, 'XTickLabel',xtlbl)   % Label Ticks
                    xlabel('Time (s)')
                    ylabel('Channel')
                    colormap(brewermap([],colormapChoice))
                    caxis(desiredPlotBoundsSubj)
                    hold on
                    for chanGoodIndSub = 1:length(notNanIndsCombined)
                        goodChanSub = notNanIndsCombined(chanGoodIndSub);
                        if ~isempty(permResultsCell{jj}.clustersSub{goodChanSub})
                            if ((~isempty(permResultsCell{jj}.clustersSub{goodChanSub}{1})) & (permResultsCell{jj}.pValuesSub{goodChanSub}(1) < 0.05))
                                rectangle('position', [tPlot(min(permResultsCell{jj}.clustersSub{goodChanSub}{1})), chanGoodIndSub-0.5,tPlot(max(permResultsCell{jj}.clustersSub{goodChanSub}{1}))-tPlot(min(permResultsCell{jj}.clustersSub{goodChanSub}{1})), 1], 'edgecolor', 'k','linewidth',2)
                            end
                        end
                    end
                    
                    currentFig = gcf;
                    currentFig.Position = [839 109 1408 1229];
                    title(['Subject ' subjName  ' Average spikes - 8 mm rereference'])
                    set(gca,'FontSize',20)
                    
                    if savePlots
                        exportgraphics(currentFig,fullfile(folderFigures,['subsample_' subjName '_perm.png']),'Resolution',600)
                        exportgraphics(currentFig,fullfile(folderFigures,['subsample_' subjName '_perm.eps']))
                    end
                    
                    figure2x2 = figure;
                    tiledlayout(2,2,'TileSpacing','Compact','Padding','Compact');
                    figure(figure2x2)
                    figure2x2.Position = [839 109 1408 1229];
                    
                    nexttile
                    imagesc([min(tPlot) max(tPlot)],[1,numPlotVec],meanSpikesReferencedGoodOnly(1:end-1,notNanIndsCombined)')
                    colormap(brewermap([],colormapChoice))
                    caxis(desiredPlotBoundsSubj)
                    hold on
                    for chanGoodInd = 1:length(notNanIndsCombined)
                        goodChan = notNanIndsCombined(chanGoodInd);
                        if ~isempty(permResultsCell{jj}.clusters{goodChan})
                            if ((~isempty(permResultsCell{jj}.clusters{goodChan}{1})) & (permResultsCell{jj}.pValues{goodChan}(1) <0.05))
                                rectangle('position', [tPlot(min(permResultsCell{jj}.clusters{goodChan}{1})), chanGoodInd-0.5,tPlot(max(permResultsCell{jj}.clusters{goodChan}{1}))-tPlot(min(permResultsCell{jj}.clusters{goodChan}{1})), 1], 'edgecolor', 'k','linewidth',2)
                            end
                        end
                    end
                    title(['Subject ' subjName ' Average spikes'])
                    set(gca,'FontSize',14)
                    
                    
                    figure(figure2x2)
                    nexttile
                    imagesc([min(tPlot) max(tPlot)],[1,numPlotVec],meanSpikesReferencedGoodOnly(1:end-1,notNanIndsCombined)')
                    colormap(brewermap([],colormapChoice))
                    caxis(desiredPlotBoundsSubj)
                    hold on
                    for chanGoodInd = 1:length(notNanIndsCombined)
                        goodChan = notNanIndsCombined(chanGoodInd);
                        if ~isempty(permResultsCell{jj}.clusters{goodChan})
                            if ((~isempty(permResultsCell{jj}.clusters{goodChan}{1})) & ( permResultsCell{jj}.adjpValuesTotal{goodChan}(1)<0.05))
                                rectangle('position', [tPlot(min(permResultsCell{jj}.clusters{goodChan}{1})), chanGoodInd-0.5,tPlot(max(permResultsCell{jj}.clusters{goodChan}{1}))- tPlot(min(permResultsCell{jj}.clusters{goodChan}{1})), 1], 'edgecolor', 'k','linewidth',2)
                            end
                        end
                    end
                    title(['Subject ' subjName ' Average spikes - FDR'])
                    set(gca,'FontSize',14)
                    
                    
                    figure(figure2x2)
                    nexttile
                    imagesc([min(tPlot) max(tPlot)],[1,numPlotVec],meanSpikesReferencedGoodOnly(1:end-1,notNanIndsCombined)')
                    colormap(brewermap([],colormapChoice))
                    caxis(desiredPlotBoundsSubj)
                    hold on
                    for chanGoodIndSub = 1:length(notNanIndsCombined)
                        goodChanSub = notNanIndsCombined(chanGoodIndSub);
                        if ~isempty(permResultsCell{jj}.clustersSub{goodChanSub})
                            if ((~isempty(permResultsCell{jj}.clustersSub{goodChanSub}{1})) & (permResultsCell{jj}.pValuesSub{goodChanSub}(1) < 0.05))
                                rectangle('position', [tPlot(min(permResultsCell{jj}.clustersSub{goodChanSub}{1})), chanGoodIndSub-0.5,tPlot(max(permResultsCell{jj}.clustersSub{goodChanSub}{1}))-tPlot(min(permResultsCell{jj}.clustersSub{goodChanSub}{1})), 1], 'edgecolor', 'k','linewidth',2)
                            end
                        end
                    end
                    title(['Subject ' subjName ' Average spikes - subsampled'])
                    set(gca,'FontSize',14)
                    
                    figure(figure2x2)
                    nexttile
                    imagesc([min(tPlot) max(tPlot)],[1,numPlotVec],meanSpikesReferencedGoodOnly(1:end-1,notNanIndsCombined)')
                    colormap(brewermap([],colormapChoice))
                    caxis(desiredPlotBoundsSubj)
                    hold on
                    for chanGoodIndSub = 1:length(notNanIndsCombined)
                        goodChanSub = notNanIndsCombined(chanGoodIndSub);
                        if ~isempty(permResultsCell{jj}.clustersSub{goodChanSub})
                            if ((~isempty(permResultsCell{jj}.clustersSub{goodChanSub}{1})) & (permResultsCell{jj}.adjpValuesTotalSub{goodChanSub}(1)<0.05))
                                rectangle('position', [tPlot(min(permResultsCell{jj}.clustersSub{goodChanSub}{1})), chanGoodIndSub-0.5,tPlot(max(permResultsCell{jj}.clustersSub{goodChanSub}{1}))-tPlot(min(permResultsCell{jj}.clustersSub{goodChanSub}{1})), 1], 'edgecolor', 'k','linewidth',2)
                            end
                        end
                    end
                    xlabel('Time (s)')
                    ylabel('Channel')
                    title(['Subject ' subjName ' Average spikes - subsampled - FDR'])
                    set(gca,'FontSize',14)
                    
                    if savePlots
                        exportgraphics(figure2x2,fullfile(folderFigures,[subjName '_2x2.png']),'Resolution',600)
                        exportgraphics(figure2x2,fullfile(folderFigures,[subjName '_2x2.eps']))
                    end
                    
                    
                end
                
            end
            
            if doConn
                disp(['Doing connected cluster permutation test'])
                
                % permutation test for background cluster based
                [clustersConn, pValuesConn, tSumsConn, permutationDistributionConn ] = permutest(referencedDataGoodOnly(1:end-1,squeeze(notNan(256,:,1)),:),referencedDataBaselineGoodOnly(:,squeeze(notNan(256,:,1)),:),false,[],[],twoSidedPerm);
                [clustersSubConn, pValuesSubConn, tSumsSubConn, permutationDistributionSubConn]  = permutest(referencedDataSubSampleGoodOnly(1:end-1,squeeze(notNanSub(256,:,1)),:),referencedDataBaselineSubSampleGoodOnly(:,squeeze(notNanSub(256,:,1)),:),false,[],[],twoSidedPerm);
                
                permResultsCellConn{jj}.subj = subjName;
                permResultsCellConn{jj}.clusters  = clustersConn;
                permResultsCellConn{jj}.pValues = pValuesConn;
                permResultsCellConn{jj}.tSums = tSumsConn;
                permResultsCellConn{jj}.permutationDistribution = permutationDistributionConn;
                permResultsCellConn{jj}.clustersSub  = clustersSubConn;
                permResultsCellConn{jj}.pValuesSub = pValuesSubConn;
                permResultsCellConn{jj}.tSumsSub = tSumsSubConn;
                permResultsCellConn{jj}.permutationDistributionSub = permutationDistributionSubConn;
                
                % extract clusters and tsums
                maxClusterConn{jj}.subj = subjName;
                maxClusterConn{jj}.clusterMax = clustersConn{1};
                maxClusterConn{jj}.tSumsMax = tSumsConn(1);
                [row,col,zdir] = ind2sub(size(referencedDataGoodOnly(1:end-1,:,:)),clustersConn{1});
                maxClusterConn{jj}.clusterMaxMatrixInds = [row col zdir];
                [boundInd] = boundary(row,col) ;
                maxClusterConn{jj}.boundInd = boundInd;
                
                maxClusterConn{jj}.clusterMaxSub = clustersSubConn{1};
                maxClusterConn{jj}.tSumsMaxSub = tSumsSubConn(1);
                [row,col,zdir] = ind2sub(size(referencedDataSubSampleGoodOnly(1:end-1,:,:)),clustersSubConn{1});
                maxClusterConn{jj}.clusterMaxMatrixIndsSub = [row col zdir];
                [boundInd] = boundary(row,col) ;
                maxClusterConn{jj}.boundIndSub = boundInd;
                disp(['Connected Cluster Permutation test done for Subject ' subjName])
                
            end
        end
        %% SVM
        if doSVM
            disp(['Doing SVM'])
            % SVM on inter vs intra
            chanListSVM = [];
            for chanGoodInd = 1:length(notNanIndsCombined)
                goodChan = notNanIndsCombined(chanGoodInd);
                if ~isempty(permResultsCell{jj}.adjpValuesTotal{goodChan})
                    if ((~isempty(permResultsCell{jj}.clusters{goodChan}{1})) & ( permResultsCell{jj}.adjpValuesTotal{goodChan}(1)<0.05))
                        chanListSVM = [chanListSVM goodChan];
                    end
                end
            end
            
            if ~isempty(chanListSVM)
                try
                    SVMModel = fitcsvm(reshape(permute(referencedDataGoodOnly(:,chanListSVM,:),[3 1 2]),length(trialsSubj),[]),trialsSubj);
                    SVMModelCrossval = crossval(SVMModel);
                    lossCrossVal = kfoldLoss(SVMModelCrossval);
                    mdlSVM = fitPosterior(SVMModel);
                    [label,scoreSvm] = resubPredict(mdlSVM);
                    if length(trialsSubj)==length(label)
                        figure
                        cm = confusionchart(trialsSubj,label);
                        %    [Xsvm,Ysvm,Tsvm,AUCsvm] = perfcurve(trialsSubj,scoreSvm(:,mdlSVM.ClassNames),'true');
                    end
                catch
                    SVMModel = nan;
                    
                end
                
                chanListSVMsub = [];
                for chanGoodInd = 1:length(notNanIndsCombined)
                    goodChanSub = notNanIndsCombined(chanGoodInd);
                    if ~isempty(permResultsCell{jj}.adjpValuesTotalSub{goodChanSub})
                        if ((~isempty(permResultsCell{jj}.clustersSub{goodChanSub}{1})) & (permResultsCell{jj}.adjpValuesTotalSub{goodChanSub}(1)<0.05))
                            chanListSVMsub = [chanListSVMsub goodChanSub];
                        end
                    end
                end
                try
                    SVMModelSubSample = fitcsvm(reshape(permute(referencedDataSubSampleGoodOnly(:,chanListSVMsub,:),[3 1 2]),length(trialsSubj),[]),trialsSubj);
                    SVMModelSubSampleCrossval = crossval(SVMModelSubSample);
                    lossCrossValSubSample = kfoldLoss(SVMModelSubSampleCrossval);
                    mdlSVMSubSample = fitPosterior(SVMModelSubSample);
                    [labelSubSample,score_svm_subsample] = resubPredict(mdlSVMSubSample);
                    if length(trialsSubj) ==length(label)
                        figure
                        cm = confusionchart(trialsSubj,labelSubSample);
                        % [Xsvm,Ysvm,Tsvm,AUCsvm] = perfcurve(trialsSubj,scoreSvmSubsample(:,mdlSVMSubSample.ClassNames),'true');
                    end
                catch
                    SVMModelSubSample = nan;
                end
                
                disp(['SVM done for Subject ' subjName])
                
                svmCell{jj}.svm = SVMModel;
                svmCell{jj}.svmSub = SVMModelSubSample;
            end
        end
        
        clearvars notNanInds notNanIndsCombined notNan notNanInds notNanSub notNanSubInds
        
    end
    
    %% captions on figures
    
    %totalFigs
    figure(figureTotal)
    xlabel('Time (s)')
    ylabel('Channel')
    sgtitle('High-density average spikes with significant channels')
    if savePlots
        exportgraphics(figureMeanSpikes,fullfile(folderFigures,'total_spikes_significant.png'),'Resolution',600)
        exportgraphics(figureMeanSpikes,fullfile(folderFigures,'total_spikes_significant.eps'))
    end
    
    figure(figureTotalSub)
    xlabel('Time (s)')
    ylabel('Channel')
    sgtitle('Sub-sampled average spikes with significant channels')
    if savePlots
        exportgraphics(figureMeanSpikes,fullfile(folderFigures,'total_spikes_significant_sub.png'),'Resolution',600)
        exportgraphics(figureMeanSpikes,fullfile(folderFigures,'total_spikes_significant_sub.eps'))
    end
    
    % mean spikes
    figure(figureMeanSpikes)
    xlabel('Time (s)')
    ylabel('Channel')
    sgtitle('Average spikes - no rereference')
    figureMeanSpikes.Position = [839 109 1408 1229];
    
    if savePlots
        exportgraphics(figureMeanSpikes,fullfile(folderFigures,'mean_spikes.png'),'Resolution',600)
        exportgraphics(figureMeanSpikes,fullfile(folderFigures,'mean_spikes.eps'))
    end
    
    figure(figureMeanRerefSubSample)
    xlabel('Time (s)')
    ylabel('Channel')
    sgtitle('Average spikes - 8 mm  rereference')
    figureMeanRerefSubSample.Position = [839 109 1408 1229];
    if savePlots
        exportgraphics(figureMeanRerefSubSample,fullfile(folderFigures,'mean_spikes_subsample_rereference.png'),'Resolution',600)
        exportgraphics(figureMeanRerefSubSample,fullfile(folderFigures,'mean_spikes_subsample_rereference.eps'))
    end
    
    figure(figureMeanReref)
    xlabel('Time (s)')
    ylabel('Channel')
    sgtitle('Average spikes - 4 mm rereference')
    figureMeanReref.Position = [839 109 1408 1229];
    if savePlots
        exportgraphics(figureMeanReref,fullfile(folderFigures,'mean_spikes_rereference.png'),'Resolution',600)
        exportgraphics(figureMeanReref,fullfile(folderFigures,'mean_spikes_rereference.eps'))
    end
    
    % mean baseline
    figure(figureBaselineMean)
    xlabel('Time (s)')
    ylabel('Channel')
    sgtitle('Average Baseline - no rereference')
    figureBaselineMean.Position = [839 109 1408 1229];
    if savePlots
        exportgraphics(figureBaselineMean,fullfile(folderFigures,'mean_baseline.png'),'Resolution',600)
        exportgraphics(figureBaselineMean,fullfile(folderFigures,'mean_baseline.eps'))
    end
    
    figure(figureBaselineMeanRerefSubSample)
    xlabel('Time (s)')
    ylabel('Channel')
    sgtitle('Average Baseline - 8 mm  rereference')
    figureBaselineMeanRerefSubSample.Position = [839 109 1408 1229];
    if savePlots
        exportgraphics(figureBaselineMeanRerefSubSample,fullfile(folderFigures,'mean_baseline_subsample.png'),'Resolution',600)
        exportgraphics(figureBaselineMeanRerefSubSample,fullfile(folderFigures,'mean_baseline_subsample.eps'))
    end
    
    figure(figureBaselineMeanReref)
    xlabel('Time (s)')
    ylabel('Channel')
    sgtitle('Average Baseline - 4 mm rereference')
    figureBaselineMeanReref.Position = [839 109 1408 1229];
    if savePlots
        exportgraphics(figureBaselineMeanReref,fullfile(folderFigures,'mean_baseline_rereference.png'),'Resolution',600)
        exportgraphics(figureBaselineMeanReref,fullfile(folderFigures,'mean_baseline_rereference.eps'))
    end
    
    % std spikes
    figure(figureStdSpikes)
    xlabel('Time (s)')
    ylabel('Channel')
    sgtitle('Std dev spikes - no rereference')
    figureStdSpikes.Position = [839 109 1408 1229];
    if savePlots
        exportgraphics(figureStdSpikes,fullfile(folderFigures,'std_spikes.png'),'Resolution',600)
        exportgraphics(figureStdSpikes,fullfile(folderFigures,'std_spikes.eps'))
    end
    
    figure(figureStdRerefSubSample)
    xlabel('Time (s)')
    ylabel('Channel')
    sgtitle('Std dev spikes - 8 mm  rereference')
    figureStdRerefSubSample.WindowState = 'fullscreen';
    if savePlots
        exportgraphics(figureStdRerefSubSample,fullfile(folderFigures,'std_spikes_subsample.png'),'Resolution',600)
        exportgraphics(figureStdRerefSubSample,fullfile(folderFigures,'std_spikes_subsample.eps'))
    end
    
    figure(figureStdReref)
    xlabel('Time (s)')
    ylabel('Channel')
    sgtitle('Std dev spikes - 4 mm rereference')
    figureStdReref.Position = [839 109 1408 1229];
    if savePlots
        exportgraphics(figureStdReref,fullfile(folderFigures,'std_spikes_rereference.png'),'Resolution',600)
        exportgraphics(figureStdReref,fullfile(folderFigures,'std_spikes_rereference.eps'))
    end
    
    % std baseline
    figure(figureBaselineStd)
    xlabel('Time (s)')
    ylabel('Channel')
    sgtitle('Std dev Baseline - no rereference')
    figureBaselineStd.Position = [839 109 1408 1229];
    if savePlots
        exportgraphics(figureBaselineStd,fullfile(folderFigures,'std_baseline.png'),'Resolution',600)
        exportgraphics(figureBaselineStd,fullfile(folderFigures,'std_baseline.eps'))
    end
    
    figure(figureBaselineStdRerefSubSample)
    xlabel('Time (s)')
    ylabel('Channel')
    sgtitle('Std dev Baseline - 8 mm  rereference')
    figureBaselineStdRerefSubSample.Position = [839 109 1408 1229];
    if savePlots
        exportgraphics(figureBaselineStdRerefSubSample,fullfile(folderFigures,'std_baseline_subsample.png'),'Resolution',600)
        exportgraphics(figureBaselineStdRerefSubSample,fullfile(folderFigures,'std_baseline_subsample.eps'))
    end
    
    figure(figureBaselineStdReref)
    xlabel('Time (s)')
    ylabel('Channel')
    sgtitle('Std dev Baseline - 4 mm rereference')
    figureBaselineStdReref.Position = [839 109 1408 1229];
    if savePlots
        exportgraphics(figureBaselineStdReref,fullfile(folderFigures,'std_baseline_rereference.png'),'Resolution',600)
        exportgraphics(figureBaselineStdReref,fullfile(folderFigures,'std_baseline_rereference.eps'))
    end
    
    save(fullfile(folderFigures,saveNameSpecific),'permResultsCell','referencedDataBaselineGoodOnly','referencedDataGoodOnly','referencedData','referencedDataBaseline','referencedDataSubSample','referencedDataBaselineSubSample','referencedDataBaselineSubSampleGoodOnly','patientsVetted','svmCell')
    clearvars permResultsCell svmCell referencedDataBaselineGoodOnly referencedDataGoodOnly referencedData referencedDataBaseline referencedDataSubSample referencedDataBaselineSubSample referencedDataBaselineSubSampleGoodOnly
    close all
    
end