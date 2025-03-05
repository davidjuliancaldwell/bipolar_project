
function fig3_call

fig = figure(3);
set(fig, 'Position', [100, 100, 1400, 900]);
set(gcf, 'Color', 'white');
pts = {'EC175', 'EC183'};
pts_ref = {'Pt. 2', 'Pt. 4'};
xldist = [0 60];
ft=[2 5 10 20 50 100 200];
ftl=cellstr(num2str(ft'));

[pt, binz, toplot, frx, binsz, Mbp_distance, ~] = fig3_EachVsAll(pts{1});

subplot1 = subplot('Position', [0.415, 0.75, 0.18, 0.18]);
elecsbrain(pt,0,[1:256],[0 0 0],'l',0,2.2,2); alpha 1;

subplot2 = subplot('Position', [0.37, 0.25, 0.25, 0.4]);
pcolorjk(binz(1:size(toplot,2)),frx,toplot); shading flat; set(gca,'ydir','normal'); ylabel('Frequency (Hz)'); xlabel('Distance (mm)'); set(gca,'fontsize',sizeoffont); %colorbar;
title({pts_ref{1} ' - ln(power), z-scored by frequency',''},'fontweight','normal')
set(gca,'yscale','log','ytick',ft,'yticklabel',ftl); xlim(xldist); clim([-1 1]*(max(abs(clim)))); colormap(gca,cmocean('curl')); %cm=cmocean('balance',100); colormap(gca,cm(:,[2 1 3]))
clim([-4.5 4.5]);

subplot3 = subplot('Position',[0.37, 0.13, 0.3, 0.06]);
histogram(make1d(Mbp_distance),[0.001 binsz:binsz:85],'facecolor',.5*[1 1 1]); set(gca,'fontsize',12); 
%xlabel('Binned bipolar distance (mm)','fontsize',sizeoffont);
ylabel('# pairs','fontweight','normal'); axis tight; grid on; cb=colorbar; set(cb,'visible','off'); xlim(xldist); set(gca,'fontsize',sizeoffont);

clear pt binz toplot frx binsz Mbp_distance cm_distance;

[pt, binz, toplot, frx, binsz, Mbp_distance, ~] = fig3_EachVsAll(pts{2});

subplot4 = subplot('Position', [0.735, 0.75, 0.18, 0.18]);
elecsbrain(pt,0,[1:256],[0 0 0],'l',0,2.2,2); alpha 1;

subplot5 = subplot('Position', [0.7, 0.25, 0.3, 0.4]);
pcolorjk(binz(1:size(toplot,2)),frx,toplot); shading flat; set(gca,'ydir','normal'); ylabel('Frequency (Hz)'); xlabel('Distance (mm)'); set(gca,'fontsize',sizeoffont); colorbar;
title({pts_ref{2} ' - ln(power), z-scored by frequency',''},'fontweight','normal')
set(gca,'yscale','log','ytick',ft,'yticklabel',ftl); xlim(xldist); clim([-1 1]*(max(abs(clim)))); colormap(gca,cmocean('curl')); %cm=cmocean('balance',100); colormap(gca,cm(:,[2 1 3]))
clim([-4.5 4.5]);

subplot6 = subplot('Position',[0.7, 0.13, 0.3, 0.06]);
histogram(make1d(Mbp_distance),[0.001 binsz:binsz:85],'facecolor',.5*[1 1 1]); set(gca,'fontsize',12); 
%xlabel('Binned bipolar distance (mm)','fontsize',sizeoffont);
%ylabel('# pairs','fontweight','normal'); 
axis tight; grid on; cb=colorbar; set(cb,'visible','off'); xlim(xldist); set(gca,'fontsize',sizeoffont);

subplot7 = subplot('Position', [0.017, 0.13, 0.3, 0.3]); 
pcolorjk(Mbp_distance); shf; hold on; 
plot([1 256; 256 256]',[1 1;1 256],'k-','linewidth',1); 
text(1, -14, '1'); text(128, -14, '128'); text(251, -14, '256');
text(263, 5, '1'); text(263, 128, '128'); text(263, 256, '256');
text(1, -33, 'Electrode #');
text(290, -1, 'Electrode #', 'Rotation', -90);
%plot([1 1 1 128 256],[1 128 256 256 256],'k.'); %plot([1 128 256 256 256]/2,[1 1 1 128 256]/2,'k.');
axis equal off; set(gca,'ydir','reverse'); 

colormap(gca,cm_distance); colorbar('NorthOutside'); 
%title('Bipolar pair distance (mm)','fontsize',sizeoffont,'fontweight','normal')
clim([0 80]);
