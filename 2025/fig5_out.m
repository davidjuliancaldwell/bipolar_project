
%% FIGURE 5
function fig5_out

fig = figure(5);
set(fig, 'Position', [100, 100, 1400, 900]);
set(gcf, 'Color', 'white');

load('stg_Devon_183_plotdata.mat');


subplot1 = subplot('Position', [0, 0.58, 0.385, 0.365]);
pt = 'EC183';
elecsbrain(pt,0,[1:256],[0 0 0],'l',0,5,2); alpha 1;
stg=getregionelecs(pt, 'stg');
elecsbrain(pt,0,[stg],[1 1 0], 0, 0, 5, 2);

xlim_all = [0 40];
xpos = [10 20 30 40 50 60];
ypos = [5 10 20 50 100 200];

subplot2 = subplot('Position', [0.409, 0.6676, 0.17, 0.187]);
pcolorjk_djc(binz(2:size(mz_zSpeechSTG,3)+1),frx,avgStimBaseSTG); shading flat; 
ylabel('Frequency (Hz)', 'FontSize',9); xlabel('Distance (mm)','FontSize',9);
set(gca,'ydir','normal');  
%text(max(xlim)+diff(xlim)/6,mean(ylim),'ln(z-score power)','fontsize',12,'rotation',90,'horizontalalignment','center')
cbar = colorbar();
cbar.Label.FontSize=8;
loc = cbar.Position;
cbar.Label.String = 'ln(z-score power)';
title({'STG Stimulus - Baseline','(z-scored by frequency)'},'fontweight','normal', 'FontSize',10)
set(gca,'yscale','log','ytick',ft,'yticklabel',ftl);
xticks([10 20 30 40 50 60]);
caxis([-maxSTGAbs,maxSTGAbs])
cmocean('balance')
xlim(xlim_all);

for i = 1:length(xpos)
    xline(xpos(i), 'Color', [0.6 0.6 0.6]);
    yline(ypos(i), 'Color', [0.6 0.6 0.6]);
end

subplot3 = subplot('Position', [0.6344, 0.6676, 0.126, 0.187]);
pcolorjk_djc(binz(2:size(mz_zNoSTSTG,3)+1),frx,avgStimBaseClustSTG); shading flat; 
set(gca,'ydir','normal'); ylabel('Frequency (Hz)','FontSize',9); xlabel('Distance (mm)','FontSize',9); 

%text(max(xlim)+diff(xlim)/4,mean(ylim),'ln(z-score power)','fontsize',12,'rotation',90,'horizontalalignment','center')
%cbar = colorbar();
%loc = cbar.Position;
%cbar.Position = [loc(1), loc(2), loc(3), loc(4)];
%cbar.Label.String = 'ln(z-score power)';
%cbar.Label.FontSize = 8.7;
title({'STG Stimulus - Baseline', ' (z-scored by frequency)', 'Significant Differences'}, ...
    'fontweight','normal', 'FontSize',10)

set(gca,'yscale','log','ytick',ft,'yticklabel',ftl);
caxis([-maxSTGAbs,maxSTGAbs])
xticks([10 20 30 40 50 60]);
cmocean('balance')
xlim(xlim_all);
for i = 1:length(xpos)
    xline(xpos(i), 'Color', [0.6 0.6 0.6]);
    yline(ypos(i), 'Color', [0.6 0.6 0.6]);
end

subplot4 = subplot('Position', [0.8376, 0.6676, 0.126, 0.187]);
plot(binz(2:size(mz_zStimSTG_gamma,1)+1),mean(mz_zStimSTG_gamma,2),'linewidth',2, 'Color', 'm');
hold on
plot(binz(2:size(mz_zNoSTSTG_gamma,1)+1),mean(mz_zNoSTSTG_gamma,2),'linewidth',2, 'Color', 'k');


for signifClust = 1:length(pValuesSpeechSTG_gamma)
    if pValuesSpeechSTG_gamma(signifClust) <= 0.05
        sigstar([binzPlotSTG(clustersSpeechSTG_gamma{signifClust}(1)),binzPlotSTG(clustersSpeechSTG_gamma{signifClust}(end))])
        
    end
end

legend({'Stimulus','Baseline'})
legend('Location', 'southeast')
title({'STG Gamma Power', 'Stimulus vs. Baseline'},'FontWeight', 'normal', 'FontSize',10)
xlim(xlim_all);
xticks([10 20 30 40 50 60]);
xlabel('Distance (mm)', 'FontSize',9)
ylabel('Averaged power (ln z-score)', 'FontSize', 9)
grid on;
set(gca, 'GridAlpha', 0.35)
ylim([-2.5 1]);

%[0.8426, 0.6876, 0.126, 0.187]

binsz=3;
subplot5 = subplot('Position', [0.8376, 0.585, 0.126, 0.038]);
histogram(make1d(MbpdistSTG),[0.001 binsz:binsz:85],'facecolor',.5*[1 1 1]); 
set(gca,'fontsize',8); %xlabel('Binned bipolar distance (mm)','fontsize',9); 
ylabel('# pairs','fontweight','normal', 'FontSize', 9); axis tight; grid on;  
xlim(xlim_all)
xticks([10 20 30 40 50 60]);
set(gca, 'GridAlpha', 0.35)

%subplot6 = subplot('Position', [0.01, 0.255, 0.159, 0.1739]);
%img2 = imread('/home/devkrish/Desktop/ec175.png');
%imshow(img2);

subplot6 = subplot('Position', [0.03, 0.185, 0.315, 0.395]);
pt = 'EC175';
elecsbrain(pt,0,[1:256],[0 0 0],'l',0,4.5,2); alpha 1;
stg=getregionelecs(pt, 'stg');
elecsbrain(pt,0,[stg],[1 1 0], 0, 0, 4.5, 2);

clear;
xlim_all = [0 40];
xpos = [10 20 30 40 50 60];
ypos = [5 10 20 50 100 200];
load('stg_Devon_175_plotdata.mat');


subplot7 = subplot('Position', [0.409, 0.32, 0.17, 0.187]);
pcolorjk_djc(binz(2:size(mz_zSpeechSTG,3)+1),frx,avgStimBaseSTG); shading flat; 
ylabel('Frequency (Hz)', 'FontSize',9); xlabel('Distance (mm)','FontSize',9);
set(gca,'ydir','normal');  
%text(max(xlim)+diff(xlim)/6,mean(ylim),'ln(z-score power)','fontsize',12,'rotation',90,'horizontalalignment','center')
cbar = colorbar();
loc = cbar.Position;
cbar.Label.String = 'ln(z-score power)';
cbar.Label.FontSize=8;
title({'STG Stimulus - Baseline','(z-scored by frequency)'},'fontweight','normal', 'FontSize',10)
set(gca,'yscale','log','ytick',ft,'yticklabel',ftl);
caxis([-maxSTGAbs,maxSTGAbs])
xticks([10 20 30 40 50 60]);
cmocean('balance')
xlim(xlim_all);

for i = 1:length(xpos)
    xline(xpos(i), 'Color', [0.6 0.6 0.6]);
    yline(ypos(i), 'Color', [0.6 0.6 0.6]);
end


subplot8 = subplot('Position', [0.6344, 0.32, 0.126, 0.187]);
pcolorjk_djc(binz(2:size(mz_zNoSTSTG,3)+1),frx,avgStimBaseClustSTG); shading flat; 
set(gca,'ydir','normal'); ylabel('Frequency (Hz)', 'FontSize',9); xlabel('Distance (mm)', 'FontSize',9); 
%text(max(xlim)+diff(xlim)/4,mean(ylim),'ln(z-score power)','fontsize',12,'rotation',90,'horizontalalignment','center')

%cbar = colorbar();
%loc = cbar.Position;
%cbar.Position = [loc(1), loc(2), loc(3), loc(4)];
%cbar.Label.FontSize=8.7;
%cbar.Label.String = 'ln(z-score power)';
title({'STG Stimulus - Baseline', '(z-scored by frequency)', 'Significant Differences'}, ...
    'fontweight','normal', 'FontSize',10)

set(gca,'yscale','log','ytick',ft,'yticklabel',ftl);
caxis([-maxSTGAbs,maxSTGAbs])
xticks([10 20 30 40 50 60]);
cmocean('balance')
xlim(xlim_all);

for i = 1:length(xpos)
    xline(xpos(i), 'Color', [0.6 0.6 0.6]);
    yline(ypos(i), 'Color', [0.6 0.6 0.6]);
end

subplot9 = subplot('Position', [0.8376, 0.32, 0.126, 0.187]);
plot(binz(2:size(mz_zStimSTG_gamma,1)+1),mean(mz_zStimSTG_gamma,2),'linewidth',2, 'Color','m')
hold on
plot(binz(2:size(mz_zNoSTSTG_gamma,1)+1),mean(mz_zNoSTSTG_gamma,2),'linewidth',2, 'Color', 'k')

for signifClust = 1:length(pValuesSpeechSTG_gamma)
    if pValuesSpeechSTG_gamma(signifClust) <= 0.05
        sigstar([binzPlotSTG(clustersSpeechSTG_gamma{signifClust}(1)),binzPlotSTG(clustersSpeechSTG_gamma{signifClust}(end))])
        
    end
end

legend({'Stimulus','Baseline'})
legend('off');
title({'STG Gamma Power', 'Stimulus vs. Baseline'},'FontWeight', 'normal', 'FontSize',10)
xlim(xlim_all);

xlabel('Distance (mm)', 'FontSize', 9)
ylabel('Averaged power (ln z-score)', 'FontSize', 9)
grid on;
set(gca, 'GridAlpha', 0.35)
xticks([10 20 30 40 50 60]);
ylim([-2.5 1]);


binsz=3;
subplot10 = subplot('Position', [0.8376, 0.2374, 0.126, 0.038]);
histogram(make1d(MbpdistSTG),[0.001 binsz:binsz:85],'facecolor',.5*[1 1 1]); 
set(gca,'fontsize',8); %xlabel('Binned bipolar distance (mm)','fontsize',9); 
ylabel('# pairs','fontweight','normal', 'FontSize', 9); axis tight; grid on;  
xlim(xlim_all)
set(gca, 'GridAlpha', 0.35)
xticks([10 20 30 40 50 60]);


%set(gca, 'XTick', []);

%title('EC183');


%% BRAIN MAPS
%{
img = imread('/home/devkrish/Desktop/ec175.png');
figure()
imshow(img)

%% LEFT PLOT

figure();
pcolorjk_djc(binz(2:size(mz_zSpeechSTG,3)+1),frx,avgStimBaseSTG); shading flat; set(gca,'ydir','normal');  set(gca,'fontsize',14);
%text(max(xlim)+diff(xlim)/6,mean(ylim),'ln(z-score power)','fontsize',12,'rotation',90,'horizontalalignment','center')
cbar = colorbar();
title('STG Stimulus - Baseline (z-scored by frequency)','fontweight','normal')
set(gca,'yscale','log','ytick',ft,'yticklabel',ftl);
% hold on
% for jj = 1:length(clusterSigSpeech)
%    plot(clustXSpeech{jj}(boundarySigSpeech{jj}),clustYSpeech{jj}(boundarySigSpeech{jj}),'k','linewidth',3)
% end
caxis([-maxSTGAbs,maxSTGAbs])
cmocean('balance')

%% MIDDLE PLOT

figure();
pcolorjk_djc(binz(2:size(mz_zNoSTSTG,3)+1),frx,avgStimBaseClustSTG); shading flat; set(gca,'ydir','normal'); ylabel('Frequency (Hz)'); xlabel('Distance (mm)'); set(gca,'fontsize',14);
%text(max(xlim)+diff(xlim)/4,mean(ylim),'ln(z-score power)','fontsize',12,'rotation',90,'horizontalalignment','center')
cbar = colorbar();
cbar.Label.String = 'ln(z-score power)';
title('STG Stim - Baseline (z-scored by frequency) Significant differences','fontweight','normal')
set(gca,'yscale','log','ytick',ft,'yticklabel',ftl);
caxis([-maxSTGAbs,maxSTGAbs])
cmocean('balance')

%% RIGHT PLOT

figure();
plot(binz(2:size(mz_zStimSTG_gamma,1)+1),mean(mz_zStimSTG_gamma,2),'linewidth',2)
hold on
plot(binz(2:size(mz_zNoSTSTG_gamma,1)+1),mean(mz_zNoSTSTG_gamma,2),'linewidth',2)

for signifClust = 1:length(pValuesSpeechSTG_gamma)
    if pValuesSpeechSTG_gamma(signifClust) <= 0.05
        sigstar([binzPlotSTG(clustersSpeechSTG_gamma{signifClust}(1)),binzPlotSTG(clustersSpeechSTG_gamma{signifClust}(end))])
        
    end
end

legend({'Stimulus','Baseline'})
title('STG Gamma Power Stimulus vs. Baseline')
xlim(xlimsTotal)

xlabel('Distance (mm)')
ylabel('Averaged power (ln z-score)')
set(gca,'fontsize',14);

%}