function loopbipolarexpedition_2025

% Code for a loop to run all patients and save display outputs
pts={'EC133','EC175','EC181','EC183','EC186','EC187','EC196','EC219','220','EC221','EC222'};
windowstocheck=250; 
nchtocheck=256; 
binsz=2;   binz=[-1 0:binsz:85];
mDiff=[];  mb_m=[];  mARb_m=[];
for p=[1 2     4 5 6 7 8     10 11]; 
    [mDiff(p,:,:),mb_m(p,:,:),mARb_m(p,:,:),binz,frx]=bipolarexpedition_EachVsAll_2023(pts{p},nchtocheck,windowstocheck,binsz);
    save('~/Desktop/bipolar_results/loop_results.mat');
end

%%
% cm=cmocean('balance',24*2);
caxdiff=[-50 150]; cmdiff=cmocean('balance',30); cmdiff(1:10,:)=[];
ft=[2 5 10 20 50 100 200]; ftl=cellstr(num2str(ft'));

figure('color','w','position',[45 249 1728 798]); 
for i=1:size(mDiff,1); subplot(3,4,i); pcolorjk(binz,frx,squeeze(mDiff(i,:,:))); cax(i,:)=clim; title(pts{i}); 
end; 
subplot(3,4,size(mDiff,1)+1); pcolorjk(binz,frx,squeeze(mean(mDiff,1))); shading flat; cax(size(mDiff,1)+1,:)=clim; title('All (mean)')
% cax=[-60 60]; %cax=[min(cax(:,1)) max(cax(:,2))]; 
colormap(cmdiff)
for i=1:12; sp(3,4,i); 
    clim(caxdiff); shading flat; set(gca,'yscale','log','ytick',ft,'yticklabel',ftl); colorbar; xlim([2.001 75]); 
end

figure('color','w'); 
pcolorjk(binz,frx,squeeze(mean(mDiff,1))); shading flat; title('All (mean)')
colormap(cmdiff)
clim(caxdiff); shading flat; set(gca,'yscale','log','ytick',ft,'yticklabel',ftl); colorbar; xlim([2.001 20]);