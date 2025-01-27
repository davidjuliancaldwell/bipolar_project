function loopbipolarexpedition

% Code for a loop to run all patients and save display outputs
pts={'EC133','EC175','EC181','EC183','EC186','EC187','EC196','EC219','220','EC221','EC222'};
mDiff=[];  mb_m=[];  mARb_m=[];
for p=[1 2     4 5 6 7 8     10 11]; 
    [mDiff(p,:,:),mb_m(p,:,:),mARb_m(p,:,:)]=bipolarexpedition_EachVsAll_2023(pts{p},256);
    save('~/Desktop/bipolar_results/loop_results');
end

%%
cm=cmocean('balance',24*2);
ft=[2 5 10 20 50 100 200]; ftl=cellstr(num2str(ft'));

figure('color','w','position',[45 249 1728 798]); 
for i=1:size(mDiff,1); subplot(3,4,i); pcolorjk(binz,frx,squeeze(mDiff(i,:,:))); cax(i,:)=clim; title(pts{i}); 
end; 
subplot(3,4,size(mDiff,1)+1); pcolorjk(binz,frx,squeeze(mean(mDiff,1))); shading flat; cax(size(mDiff,1)+1,:)=clim; title('All (mean)')
cax=[-60 60]; %cax=[min(cax(:,1)) max(cax(:,2))]; 
colormap(cm)
for i=1:12; sp(3,4,i); clim(cax); shading flat; set(gca,'yscale','log','ytick',ft,'yticklabel',ftl); colorbar; xlim([2.001 75]); end

figure('color','w'); 
pcolorjk(binz,frx,squeeze(mean(mDiff,1))); shading flat; title('All (mean)')
colormap(cm)
clim(cax); shading flat; set(gca,'yscale','log','ytick',ft,'yticklabel',ftl); colorbar; xlim([2.001 20]);