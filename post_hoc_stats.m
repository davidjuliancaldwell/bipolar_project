% script to do statistics post hoc on data
%%
savePlots = 1;

folderDataBase = '/Users/davidcaldwell/Box/KLEENLAB/David/Results/2025/';
folderFiguresCell = {fullfile(folderDataBase,'LL20'),fullfile(folderDataBase,'LL40'),fullfile(folderDataBase,'LL100'),fullfile(folderDataBase,'absDer')};
saveName = {'LL20','LL40','LL100','absDer'};
statsStruct = struct;

outputTable = table;

% load in data
for jj = 1:length(saveName)
    processedInt = saveName{jj};
    folderFigures = folderFiguresCell{jj};
    dataFile  = fullfile(folderFigures,[processedInt '.mat']);
    load(dataFile);
    
    %%
    subjCell = [];
    for jjj = 1:length(permResultsCell)
        subjCell{end+1} = permResultsCell{jjj}.subj;
    end
    
    numEls = length(permResultsCell);
   % c= cmocean('thermal',numEls);
   c = brewermap(numEls,'Set3');
    disp(numEls);
    
    %%
    numChannelsVec = [];
    meanWidthVec = [];
    meanWidthSID = {};
    numChannelsVecSub = [];
    meanWidthVecSub = [];
    numChannelsSID = {};
    %
    for jjj = 1:length(permResultsCell)
        
        % if sum(strcmp(patientsVetted{jj},patientsInt))
        width = permResultsCell{jjj}.meanWidth;
        widthSub = permResultsCell{jjj}.meanWidthSub;
        num = permResultsCell{jjj}.numSigChannels;
        numSub = permResultsCell{jjj}.numSigChannelsSub;
        subj = permResultsCell{jjj}.subj;
        
        if ~isnan(width) & ~isnan(widthSub)
            meanWidthVec = [meanWidthVec width];
            meanWidthVecSub = [meanWidthVecSub widthSub];
            meanWidthSID{end+1} = subj;
        end
        
        if ~isnan(num) & ~isnan(numSub)
            numChannelsVec = [numChannelsVec num];
            numChannelsVecSub = [numChannelsVecSub numSub];
            numChannelsSID{end+1} = subj;
        end
        %  end
        
    end
    %% plot
    
    colorHist = cmocean('thermal',2);
    colorParallel = brewermap(numEls,'Set3');
    
    histFig = figure;
    tiledlayout(2,2,'TileSpacing','Compact','Padding','Compact');
    histFig.Position = [839 109 1408 1229];
    nexttile
    histogram(meanWidthVec,BinWidth=10)
    hold on
    histogram(meanWidthVecSub,BinWidth=10)
    legend({'High density','Sub-sampled'})
    title([processedInt ' Mean Width Histogram'])
    xlim([0 max([meanWidthVec(:);meanWidthVecSub(:)])+10])
    set(gca,'FontSize',14)
    
    nexttile
    histogram(numChannelsVec,BinWidth=1)
    hold on
    histogram(numChannelsVecSub,BinWidth=1)
    legend({'High density','Sub-sampled'})
    title([processedInt ' Number of Channels Histogram'])
    set(gca,'FontSize',14)
    xlim([0 max([numChannelsVec(:);numChannelsVecSub(:)])+3])
    
    colormapSpecificWidth = [];
    for jjj = 1:length(subjSpecific)
        ind = find(strcmp(subjCell,subjSpecific(jjj)));
        colormapSpecificWidth = [colormapSpecificWidth; colorParallel(ind,:)];
    end
    
    ax= nexttile;
    grid(ax,'on')
    xlim([0.5 2.5])
    ylim([0 max([meanWidthVec(:);meanWidthVecSub(:)])+10])
    hold on
    p1=plot([1,2],[meanWidthVec;meanWidthVecSub]','-o','LineWidth',4);
    
    xlabel('High density vs. Sub-sampled')
    ylabel('Mean width')
    title([processedInt ' mean spike subject wise comparison'])
    set(gca,'FontSize',14)
    ax.XTickLabel ={'','High density','','Sub-sampled',''};
    colororder(gca,colormapSpecificWidth)
    for colorInd = 1:length(meanWidthVec)
        p1(colorInd).MarkerFaceColor = colormapSpecificWidth(colorInd,:);
    end
    
    colormapSpecificNum = [];
    for jjj = 1:length(subjSpecific)
        ind = find(strcmp(subjCell,subjSpecific(jjj)));
        colormapSpecificNum = [colormapSpecificNum; colorParallel(ind,:)];
    end

    
    ax = nexttile;
    xlim([0.5 2.5])
    grid(ax,'on')
    ylim([0 max([numChannelsVec(:);numChannelsVecSub(:)])+3])
    hold on
    p2=plot([1,2],[numChannelsVec;numChannelsVecSub]','-o','LineWidth',4);
    ylabel('Number of channels with significant spikes')
    xlabel('High density vs. Sub-sampled')
    title([processedInt ' number spike channels subject wise comparison'])
    set(gca,'FontSize',14)
    ax.XTickLabel ={'','High density','','Sub-sampled',''};
    colororder(gca,colormapSpecificNum)
    for colorInd = 1:length(numChannelsVec)
        p2(colorInd).MarkerFaceColor = colormapSpecificNum(colorInd,:);
    end
    
    if savePlots
        exportgraphics(histFig,fullfile(folderFigures,[processedInt '_2x2.png']),'Resolution',600)
        exportgraphics(histFig,fullfile(folderFigures,[processedInt '_2x2.eps']))
    end
    
    
    %% stats
    [pNum,hNum,statsNum] = signrank(numChannelsVec,numChannelsVecSub);
    
    [pWidth,hWidth,statsWidth] = signrank(meanWidthVec,meanWidthVecSub);
    
    statsStruct.name{jj} = processedInt;
    statsStruct.meanWidthSID{jj} = meanWidthSID;
    statsStruct.numChannelsSID{jj} = numChannelsSID;
    statsStruct.pNum{jj} = pNum;
    statsStruct.pWidth{jj} = pWidth;
    statsStruct.meanWidthVec{jj} = meanWidthVec;
    statsStruct.meanWidthVecSub{jj} = meanWidthVecSub;
    statsStruct.numChannelsVec{jj} = numChannelsVec;
    statsStruct.numChannelsVecSubj{jj} = numChannelsVecSub;
    statsStruct.meanWidthDiff{jj} = meanWidthVec - meanWidthVecSub;
    statsStruct.numChannelsDiff{jj} = numChannelsVec - numChannelsVecSub;
    
end
%%
figure
tiledlayout(2,1,'TileSpacing','Compact','Padding','Compact');
grid(ax,'on')
nexttile
hold on
for jj = 1:length(saveName)
    subjSpecific = statsStruct.meanWidthSID{jj};
    colormapSpecific = [];
    for jjj = 1:length(subjSpecific)
        ind = find(strcmp(subjCell,subjSpecific(jjj)));
        colormapSpecific = [colormapSpecific; c(ind,:)];
    end
    swarmchart(jj,statsStruct.meanWidthDiff{jj},[],colormapSpecific,'filled');
end
title('Differences in Mean Spike Duration by Condition')
xticks([0 1 2 3 4 5])
xticklabels({'','LL20','LL40','LL100','Absolute Derivative',''})

nexttile
grid(ax,'on')
hold on
for jj = 1:length(saveName)
    subjSpecific = statsStruct.numChannelsSID{jj};
    colormapSpecific = [];
    for jjj = 1:length(subjSpecific)
        ind = find(strcmp(subjCell,subjSpecific(jjj)));
        colormapSpecific = [colormapSpecific; c(ind,:)];
    end
    swarmchart(jj,statsStruct.numChannelsDiff{jj},[],colormapSpecific,'filled');
end
title('Differences in Number of Channels with Spikes by Condition')
xticks([0 1 2 3 4 5])
xticklabels({'','LL20','LL40','LL100','Absolute Derivative',''})
xlabel('condition')

tempFig = gcf;
if savePlots
    exportgraphics(tempFig,fullfile(folderDataBase,['swarm_across_conditions.png']),'Resolution',600)
    exportgraphics(tempFig,fullfile(folderDataBase,['swarm_across_conditions.eps']))
end

%%

statsCell = {};
figure
tiledlayout(2,1,'TileSpacing','Compact');
ax = nexttile;
grid(ax,'on')

hold on

for jj = 1:length(saveName)
    subjSpecific = statsStruct.meanWidthSID{jj};
  
    colormapSpecific = [];
    for jjj = 1:length(subjSpecific)
        ind = find(strcmp(subjCell,subjSpecific(jjj)));
        colormapSpecific = [colormapSpecific; c(ind,:)];
        statsCell{ind}.cmap = c(ind,:);
        if jj == 1 
            statsCell{ind}.data{1} = statsStruct.meanWidthDiff{jj}(jjj);
            statsCell{ind}.cond{1} = jj;
        else
            statsCell{ind}.data{end+1} = statsStruct.meanWidthDiff{jj}(jjj);
            statsCell{ind}.cond{end+1} = jj;
                        
        end
    end
end

for jjj =  1:length(subjCell)
    if ~isempty(statsCell{jjj})
        linePlotPostHoc(jjj) = plot(cell2mat(statsCell{jjj}.cond),cell2mat(statsCell{jjj}.data),'-o','linewidth',4);
    end
end
colororder(colormapSpecific);
for colorInd = 1:length(subjCell)
    if ~isempty(statsCell{colorInd})
        linePlotPostHoc(colorInd).MarkerFaceColor = c(colorInd,:);
    end
end

title('Differences in Mean Spike Duration by Condition')
ylabel('Difference in Mean Spike Duration (ms)')
xticks([0 1 2 3 4 5])
xlim([0 5])
xticklabels({'','','','','',''})
%
ax = nexttile;
grid(ax,'on')
hold on
for jj = 1:length(saveName)
    subjSpecific = statsStruct.numChannelsSID{jj};
  
    colormapSpecificN = [];
    for jjj = 1:length(subjSpecific)
        ind = find(strcmp(subjCell,subjSpecific(jjj)));
        colormapSpecificN = [colormapSpecificN; c(ind,:)];
        statsCell{ind}.cmapN = c(ind,:);
        if jj == 1 
            statsCell{ind}.dataN{1} = statsStruct.numChannelsDiff{jj}(jjj);
            statsCell{ind}.condN{1} = jj;
        else
            statsCell{ind}.dataN{end+1} = statsStruct.numChannelsDiff{jj}(jjj);
            statsCell{ind}.condN{end+1} = jj;
                        
        end
    end
end

for jjj =  1:length(subjCell)
    if ~isempty(statsCell{jjj})
        linePlotPostHocN(jjj) = plot(cell2mat(statsCell{jjj}.condN),cell2mat(statsCell{jjj}.dataN),'-o','linewidth',4);
    end
end
colororder(colormapSpecificN);
for colorInd = 1:length(subjCell)
    if ~isempty(statsCell{colorInd})
        linePlotPostHocN(colorInd).MarkerFaceColor = c(colorInd,:);
    end
end



title('Differences in Number of Spikes Detected by Condition')
ylabel('Number of channels')
xlabel('Condition')
xticks([0 1 2 3 4 5])
xlim([0 5])
xticklabels({'','LL20','LL40','LL100','Absolute Derivative',''})

tempFig = gcf;
tempFig.Position = [986 636 581 702];
if savePlots
    exportgraphics(tempFig,fullfile(folderDataBase,['across_conditions.png']),'Resolution',600)
    exportgraphics(tempFig,fullfile(folderDataBase,['across_conditions.eps']))
end
%


%%
saveNameSpecific = 'spikeStats.mat';
save(fullfile(folderDataBase,saveNameSpecific),'statsStruct');




