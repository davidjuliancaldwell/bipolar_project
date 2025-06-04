
function loopbipolarexpedition(g1d2s3)

% g1d2s3=3;


none1sqrt2log3=2; % 1: no transform, 2: square root, 3: log

%function [mDiff, mb_m, mARb_m, binz, frx] = d_bipolarexpedition_EachVsAll_2023(pt, nchtocheck, windowstocheck)

% Code for a loop to run all patients and save display outputs
if g1d2s3==1 %grids
    pts={'EC133','EC175','EC183','EC186','EC187','EC196','EC219','EC221','EC222'}; %have grids (181, 220 don't)
    component_type='grids';
elseif g1d2s3==2 %depths
    pts={'EC133','EC175','EC181','EC183','EC186','EC187','EC196','EC219','EC220','EC221','EC222'}; %all have depths
    component_type='depths';
elseif g1d2s3==3 %strips
    pts={'EC133','EC175','EC183','EC186','EC187','EC196','EC219','EC221','EC222'}; %have strips (181, 220 don't)
    component_type='strips';
end
%depth_check = [320, 340, 84, 298, 318, 318, 308, 468, 84, 404, 324]; %depth elec nums
mDiff=[];  mb_m=[];  mARb_m=[]; Mbp_distance = {}; %rec_lens = {};

for p = 1:size(pts, 2)
    %rec_lens{end+1} = devon_EachVsAll_cleaned(pts{p});
    [mDiff(p,:,:),mb_m(p,:,:),mARb_m(p,:,:), Mbp_distance{p}]= EachVsAll_cleaned_June3(pts{p},none1sqrt2log3,g1d2s3); %concat, add extra dim
    disp(['FINISHED LOADING: ' pts{p}]);
    saveas(gcf,['/Users/jonathankleen/Desktop/bipolar_results/' component_type '_June3_2025_3mm/' pts{p} '___ex.png']); close
    tic; save(['/Users/jonathankleen/Desktop/bipolar_results/' component_type '_June3_2025_3mm/temp_data_loop.mat']); toc
end

%%

mDiff_re = mDiff;
mb_m_re = mb_m;
mARb_m_re = mARb_m;
distance = Mbp_distance;

save(['/Users/jonathankleen/Desktop/bipolar_results/' component_type '_June3_2025_3mm/m_EachVsAll_' component_type '.mat'], 'mDiff_re', 'mb_m_re', 'mARb_m_re');
save(['/Users/jonathankleen/Desktop/bipolar_results/' component_type '_June3_2025_3mm/' component_type '_distances.mat'], 'distance');

end

%save('depths_etc.mat', 'binz_re', 'frx_re');

%%

%{
cm=cmocean('balance',24*2);
ft=[2 5 10 20 50 100 200]; ftl=cellstr(num2str(ft'));

figure('color','w','position',[45 249 1728 798]); 
for i=1:size(mDiff,1); subplot(3,4,i); pcolorjk(binz,frx,squeeze(mDiff(i,:,:))); cax(i,:)=clim; title(pts{i}); 
end; 
subplot(3,4,size(mDiff,1)+1); pcolorjk(binz,frx,squeeze(mean(mDiff,1))); shading flat; cax(size(mDiff,1)+1,:)=clim; title('All (mean)')
cax=[-60 60]; %cax=[min(cax(:,1)) max(cax(:,2))]; 
colormap(cm)
for i=1:12; sp(3,4,i); clim(cax); shading flat; set(gca,'yscale','log','ytick',ft,'yticklabel',ftl); colorbar; xlim([2.001 52]); end

figure('color','w'); 
pcolorjk(binz,frx,squeeze(nanmean(mDiff,1))); shading flat; title('All (mean)')
colormap(cm)
clim(cax); shading flat; set(gca,'yscale','log','ytick',ft,'yticklabel',ftl); colorbar; xlim([4.001 20]);

%}