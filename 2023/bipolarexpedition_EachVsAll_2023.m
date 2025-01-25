function [mDiff,mb_m,mARb_m,binz,frx]=bipolarexpedition_EachVsAll_2023(pt,nchtocheck,windowstocheck)
% BIPOLAR PAIR ANALYSIS: EACH VS. ALL
% see loopbipolarexpedition.m to loop across patients and analyze
% EC175 and EC183 both have intact 16x16 square grids (channel #s 1:256)

if ~exist('pt','var')||isempty(pt); pt='EC175'; end %pt='EC175'; % EC175 and EC183 both have intact 16x16 square grids (channel #s 1:256)
if ~exist('nchtocheck','var')||isempty(nchtocheck); nchtocheck=128*2; end
if ~exist('windowstocheck','var')||isempty(windowstocheck); windowstocheck=250; end %each window is 1 second of data (non-overlapping)

none1sqrt2log3=3; % 1: no transform, 2: square root, 3: log
g1s2d3=1; % use either grids (1) or strips (2) or depths (3) but not the others
binsz=3; % bin size in mm

xldist=[0 70];
doanglerange=0;
sizeoffont=12;




cm=cool(6); cm(1,:)=[0 0 0];
datadir=['/Volumes/KLEEN_DRIVE/bipolar_expedition/'];
ptdatadir=[datadir 'baseline-high-density-data/'];
u=dir(ptdatadir); uptbl={}; for i=1:length(u); uname=u(i).name; uptbl{i,1}=uname(1:end-28); end; uptbl(1:2)=[]; clear i u uname

load([datadir '/taggedspikes_April2022']);
sfx=512;
frxrange=[2 200]; %frequency range to examine
ft=[2 5 10 20 50 100 200]; ftl=cellstr(num2str(ft')); %frequency labels for plots

p=find(strcmpi(pts,pt)); %patient number ID
pblocks=strfind(uptbl,pts{p}); 
for i=1:length(pblocks); 
    isbl(i,1)=~isempty(pblocks{i}); 
end
ptbl=find(isbl); if ~isempty(ptbl); disp(['Loading ' pts{p} ' blocks...']); end





% load all blocks for this patient and stack their baseline windows together
d=[]; nwind=0;
for b=1:length(ptbl); disp(uptbl{ptbl(b)})
    % load using using "_jk" versions of baseline windows, updated 2/2022
    ptpath = fullfile(ptdatadir, [uptbl{ptbl(b)} '_baselineWindows_fromraw.mat']);
    load(ptpath);
    % get rid of baseline windows containing spikes or artifact
    spksarti=hasspk | hasarti;
    nonspks_windows(:,spksarti)=[];
    hasstim(spksarti)=[]; %update indices for which windows overlap with stimuli/speech
    hasspeech(spksarti)=[]; 


    clear hasspkvec hasspk hasartivec hasarti spksarti % now clear spike- and artifact-related variables from workspace

    % convert to 3D matrix, combine all windows from consecutive blocks for each patient
    for i=1:size(nonspks_windows,2)
        d(:,:,i+nwind)=nonspks_windows{2,i}'; 
    end
    nwind=size(d,3);

    clear nonspks_windows info
end; clear b




nch=size(d,2); 

% load electrode component infor (grid/strip/depth and how many linear contacts they have in a row
% [bpN,bpT]=xlsread(['/Users/davidcaldwell/code/high_density_ecog/AN_ElectrodeInfoTDT.xlsx'],pts{p});

an_electrode_info_path = fullfile(data_root, 'bipolar_expedition', 'AN_ElectrodeInfoTDT.xlsx');
[bpN,bpT]=xlsread(an_electrode_info_path, pts{p});

[em,eleclabels,anatomy]=getelecs(pts{p},2);













%% look at either grids, strips, or depths, and nan the others
          for r=1:size(bpT,1)
            if [g1s2d3~=1 && any(strcmpi(bpT(r,2),{'grid','minigrid'}))]  || ...
               [g1s2d3~=2 &&     strcmpi(bpT(r,2),'strip')]               || ...
               [g1s2d3~=3 &&     strcmpi(bpT(r,2),'depth')];
                       d(:,bpN(r,1):bpN(r,2),:)=nan; %nan out component that isn't relevant to this run (see g1s2d3 above)
              bp_distance(bpN(r,1):bpN(r,2))  =nan; %nan out their corresponding distances (irrelevant for this run)
              bp_angle   (bpN(r,1):bpN(r,2))  =nan; %and angles, similarly
            end
          end; clear r


%% bad channels
badchI=isnan(mean(mean(d,1),3))'; %any channel with nans in any window will be a "bad channel"
badchI(badchidx{p})=1; %follow up with all marked as bad channels in original preprocessing
x=[]; xbch=true(size(badchI)); for i=1:size(bpT,1); x=[x bpN(i,1):bpN(i,2)]; end; xbch(x)=false; 
badchI(xbch)=1; %find and nan empty channels unaccounted for by component rows
d(:,badchI,:)=nan;
okc=~badchI; clear x xbch


%*at this point, isolate speech or stim or non-speech/stim

% %*at this point, isolate STG or IFG channels
% [STG]=getelecs_region(pt,'stg',2);
% [IFG]=getelecs_region(pt,{'po','pt'},2);

windowstocheck=min([windowstocheck size(d,3)]);
windowstocheck=1:windowstocheck; %convert to a vector of windows, 1:X


%% ALL PAIRS (each vs. all others) analysis and example plot
 d=d(:,:,windowstocheck); clear Straces_allch; %free up RAM by getting rid of whatever won't be used (only using first ___ number of windows)
                 % ***opportunity here to select speech or stim windows
 [M,M_averef,Mbp_distance,frx,~,Mbp_angle]=bpspectra_EachVsAll_2023(d,sfx,frxrange,em,nchtocheck,none1sqrt2log3);
     % M is bipolarchannel1 X bipolarchannel2 X frx X 1secwindow
 nfrx=length(frx);
  
%  % average across frequencies to create log scale bins
%  frxbins=[2,3,4,5,6,7,8,10,12,15,18,22,27,33,41,50,60,74,90,110,134,163,198];
 
%  %% mean across windows
%  M=squeeze(mean(M,4)); 
%  MRA=squeeze(mean(MRA,4)); 


%% ANGLE RANGE: if desired, subselect bipolar angle within a given range
if doanglerange
    m=M; Md=Mbp_distance; Ma=Mbp_angle; %make a copy (only do this once) in case you want to repeat later with different angle range (use next chunk of code)
    
    anglemin= 45; anglemax= 90;
    M=m; Mbp_distance=Md; Mbp_angle=Ma; 
    ww=0;
    for c1=1:size(M,1);
    for c2=1:size(M,2);
      if rad2ang(Mbp_angle(c1,c2))<anglemin || rad2ang(Mbp_angle(c1,c2))>anglemax;
        Mbp_distance(c1,c2,:,:)=nan;
        Mbp_angle(c1,c2,:,:)=nan;
        disp([num2str(c1) ' ' num2str(c2)])
        ww=ww+1;
      end
    end
    end
    ww
end

%% Unstack and line up into 2D matrix to create channels^2 X frequencies for easier indexing-->binning
mb=[]; 
mARb=[]; 
binz=[-1 0:binsz:85]; % can change binsz PRN at this point to test different bin resolutions on the plots below
clear mx my mz  
nbinz=length(binz)-1;
 parfor w=windowstocheck; disp(num2str(w)); % parfor here to run the windows
    Mflat=[];
      Mbp_bpflat_distance=[];
      M_bpflat_angle=[];
    M_averef_flat=[];
    for c2=1:nchtocheck; 
        Mflat=[Mflat; squeeze(M(:,c2,:,w))];               
        Mbp_bpflat_distance=[Mbp_bpflat_distance; Mbp_distance(:,c2)]; % corresponding distance index
        M_bpflat_angle   =[M_bpflat_angle;    Mbp_angle(:,c2)]; % corresponding angle index
        M_averef_flat=[M_averef_flat; squeeze(M_averef(:,c2,:,w))];              
    end
    
    %bin by distance and take the mean, creating: frequency X binned distance
       % first for referential (distance = 0)
    for i=1:nbinz % binz will be including >lower bound and up to and including (<=) upper bound
        mb  (:,i,w)=nanmean(Mflat        (Mbp_bpflat_distance>binz(i) & Mbp_bpflat_distance<=binz(i+1),:),1); 
        mARb(:,i,w)=nanmean(M_averef_flat(Mbp_bpflat_distance>binz(i) & Mbp_bpflat_distance<=binz(i+1),:),1); %AR means Average Reference
    end
    
    disp([num2str(round(windowstocheck(w)/windowstocheck(end)*100,1)) '% of windows'])
end; disp('Done')

binz(1)=[]; %remove negative 1, only there to get referential channels (distance=0) for the parfor loop

%         mb(:,1,w)  =nanmean(Mflat      (Mflatbp_distance==0,:),1);
%         mARb(:,1,w)=nanmean(Mrefaveflat(Mflatbp_distance==0,:),1);




%% --- Figure 3 --- (Omnidirectional re-referencing)

% % NOTES:
% % mb is frequency X binned distance X windows, and you can average across windows
% % binindex_min_ltmax tells you for each column of mb what was the 
% %         1] minimum (>0) distance, and
% %         2] the "less than max" (<) distance
% %         that were used to index bipolar pairs for that bin
% %         Note: first bin includes zero, which corresponds to
% %         the bin containing the referential channels


% % zscore the log transformed power according to frequency

% transform power before calculating 
switch none1sqrt2log3
  case 2
        mb__m=squeeze(mean(sqrt(mb),3)); 
        mARb__m=squeeze(mean(sqrt(mARb),3)); 
         nfrx=length(frx);
         mb__m_z=mb__m; mARb__m_z=mARb__m;
         for i=1:nfrx; nns=~isnan(mb__m(i,:)); 
             mb__m_z(i,nns)=zscore((mb__m(i,nns))); 
             mARb__m_z(i,nns)=zscore((mARb__m(i,nns))); 
         end; 
  case 3
        mblogm=squeeze(mean(log(mb),3)); 
        mARblogm=squeeze(mean(log(mARb),3)); 
         nfrx=length(frx);
         mb__m_z=mblogm; 
         mARb__m_z=mARblogm;
         for i=1:nfrx; nns=~isnan(mblogm(i,:)); 
             mb__m_z(i,nns)=zscore((mblogm(i,nns))); 
             mARb__m_z(i,nns)=zscore((mARblogm(i,nns))); 
         end; 
end


 figure('color','w','position',[[54 223 1907 1102]]); colormap(parula); 
 
 % plot log-transformed version
 toplot=mean((mb__m_z),3);  % mb or mb_z
 cm_distance=flipud(cmocean('deep',1+ceil(max(max(Mbp_distance)))));

% Plot confusion matrix of bipolar distances between all pairs
 subplot(8,3,1:3:7); pcolorjk(Mbp_distance); shf; hold on; %plot([1 1; 1 256]',[1 256;256 256]','k-','linewidth',1); plot([1 1 1 128 256],[1 128 256 256 256],'k.'); %plot([1 128 256 256 256]/2,[1 1 1 128 256]/2,'k.'); 
    axis equal off; set(gca,'ydir','normal','xdir','reverse'); colormap(gca,cm_distance); colorbar('fontsize',12); title('Bipolar pair distance (mm)','fontsize',sizeoffont,'fontweight','normal')
    clim([0 80]); 
% Plot confusion matrix of bipolar distances between all pairs
 subplot(8,3,2:3:8); pcolorjk(rad2deg(Mbp_angle)); shf; hold on; %plot([1 1; 1 256]',[1 256;256 256]','k-','linewidth',1); plot([1 1 1 128 256],[1 128 256 256 256],'k.'); %plot([1 128 256 256 256]/2,[1 1 1 128 256]/2,'k.'); 
    axis equal off; set(gca,'ydir','normal','xdir','reverse'); colormap(gca,hsv(360)); clim([-180 180]); colorbar('fontsize',12); title('Bipolar pair angle (deg)','fontsize',sizeoffont,'fontweight','normal')

% Histogram of number of pairs per bin
 subplot(8,2,7); histogram(make1d(Mbp_distance),0:binsz:85,'facecolor',.5*[1 1 1]); set(gca,'fontsize',12); xlabel('Binned bipolar distance (mm)','fontsize',sizeoffont); 
 ylabel('Counts/bin','fontweight','normal'); axis tight; grid on; cb=colorbar; set(cb,'visible','off'); xlim(xldist)


% Histogram of angles
if doanglerange
    subplot(8,6,22); angs=make1d(Mbp_angle); nns=~isnan(angs);
 polarhistogram(angs(nns),-pi:2*pi*(5/360):pi,'facecolor',.5*[1 1 1]); title('Bipolar angle distribution (5^o steps)','fontsize',sizeoffont,'fontweight','normal'); 
 maxrt=max(get(gca,'rtick')); 
 set(gca,'rtick',[min(get(gca,'rtick')) maxrt/2 maxrt],'ThetaTick',[0 90 180 270],'fontSize',9); 
end

if doanglerange;  
% plot illustration of location of bipolar pairs, colored by same distance scale
 subplot(2,3,3); hold on; 
 elecsbrain(pt,0,[1:nchtocheck],[0 0 0],'l',0,5,2); alpha 0.05; litebrain('r',0); zoom(1.5)
 for c1=1:size(Mbp_angle,1)
   for c2=1:size(Mbp_angle,2)
       if ~isnan(Mbp_angle(c1,c2))
            plot3([em(c1,1) em(c2,1)],[em(c1,2) em(c2,2)],[em(c1,3) em(c2,3)],'-','color',cm_distance(1+round(Mbp_distance(c1,c2)),:),'LineWidth',.5)%,'LineWidth',.75)
       end
   end
 end
end

 subplot(2,2,3); 
 pcolorjk(binz(1:size(toplot,2)),frx,toplot); shading flat; set(gca,'ydir','normal'); ylabel('Frequency (Hz)'); xlabel('Distance (mm)'); set(gca,'fontsize',sizeoffont); colorbar; 
 title({pt ' - ln(power), z-scored by frequency',''},'fontweight','normal')
 set(gca,'yscale','log','ytick',ft,'yticklabel',ftl); xlim(xldist); clim([-1 1]*(max(abs(clim)))); colormap(gca,cmocean('curl')); %cm=cmocean('balance',100); colormap(gca,cm(:,[2 1 3]))
 clim([-4.5 4.5]);

 subplot(2,2,4)
 [mx,my]=meshgrid(binz(1:size(toplot,2))+binsz/2,frx);
 surf(mx,my,(toplot)); xlabel('Bipolar distance (mm)','fontsize',sizeoffont); ylabel('Frequency (Hz)','fontsize',sizeoffont); zlabel({'ln(power),','z-scored by frequency'},'fontsize',sizeoffont)
 view(140,25); set(gca,'xdir','reverse','ydir','reverse','ytick',ft,'yticklabel',ftl,'fontsize',sizeoffont); 
 colormap(gca,cmocean('curl')); clim([-1 1]*(max(abs(clim)))); 
 xlim([binsz/2 85]);  axis tight; 
 set(gca,'yscale','log'); 
 xlim(xldist)
 ylim([0 200]); 
 zlim([-4.5 4.5]);
 if doanglerange; 
     title([num2str(anglemin) ' to ' num2str(anglemax)]); 
 end


if doanglerange %loop through consecutive angle ranges of bipolar pair orientations to see effects
 loopangles
 return
end

%% bipolar power minus mean referential power for all pairs
%% --- Figure 4 --- (relative changes in power)
% will also add a line at 10mm for visualization of this common clinical inter-electrode distance

% transform power before calculating perecent change
mb_=mb; mARb_=mARb; % make a copy... then transform the copy
switch none1sqrt2log3
    case 1
        txtyp='raw'; % no transform, power in raw form
        mb_(~isnan(mb_))      =    (mb_(~isnan(mb_)));
        mARb_(~isnan(mARb_))  =    (mARb_(~isnan(mARb_)));
    case 2
        txtyp='square root'; % square root transform
        mb_(~isnan(mb_))      =sqrt(mb_(~isnan(mb_)));
        mARb_(~isnan(mARb_))  =sqrt(mARb_(~isnan(mARb_)));
    case 3
        txtyp='natural log'; % log transform --> problematic because taking % change of certain small negative values creates extreme values
        mb_(~isnan(mb_))      =log (mb_(~isnan(mb_)));
        mARb_(~isnan(mARb_))  =log (mARb_(~isnan(mARb_)));
end

%     % zscore acrosss freqs and windows for each distance bin
%     for i=1:nbinz; I=mb_  (:,i,:); I(~isnan(I))=zscore(I(~isnan(I))); mb_  (:,i,:)=I; end
%     for i=1:nbinz; I=mARb_(:,i,:); I(~isnan(I))=zscore(I(~isnan(I))); mARb_(:,i,:)=I; end
    
% mean across windows
mb_m=squeeze(mean(mb_,3)); 
mARb_m=squeeze(mean(mARb_,3)); 

psig=0.05;

xl=[binsz*2 binz(end)]; fl=frx([1 end]);

figure('color','w','position',[315 1 1486 1046]); colormap(jet); clear cax; %
 subplot(3,2,1)
 pcolorjk(binz(1:size(mb_m,2))-binsz,frx,mARb_m); shading flat; ylabel('Frequency (Hz)'); xlabel('Distance (mm)'); set(gca,'fontsize',sizeoffont); colorbar; 
 text(max(xlim)+diff(xlim)/4,mean(ylim),'(power)','fontsize',12,'rotation',90,'horizontalalignment','center')
 title({'Referential',['(average of pair)' ', ' txtyp]},'fontweight','normal')
 set(gca,'ydir','normal','yscale','log','ytick',ft,'yticklabel',ftl,'xlim',xl,'ylim',fl);
 cax(1,:)=clim;

 subplot(3,2,2)
 pcolorjk(binz(1:size(mb_m,2))-binsz,frx,mb_m); shading flat; ylabel('Frequency (Hz)'); xlabel('Distance (mm)'); set(gca,'fontsize',sizeoffont); colorbar; 
 text(max(xlim)+diff(xlim)/4,mean(ylim),'(power)','fontsize',12,'rotation',90,'horizontalalignment','center')
 title(['Bipolar' ', ' txtyp],'fontweight','normal')
 set(gca,'ydir','normal','yscale','log','ytick',ft,'yticklabel',ftl,'xlim',xl,'ylim',fl);
 cax(2,:)=clim;
 cax=(([min(cax(:,1),[],1) max(cax(:,2),[],1)])); % get extreme min and max of the two conditions
 subplot(3,2,1);   clim(cax); xlim(xldist); yline(10,'k-',.75); 
 subplot(3,2,2);   clim(cax); xlim(xldist); yline(10,'k-',.75); 
 
 subplot(3,2,3)
 mDiff=((mb_m-mARb_m)./mARb_m)*100; difftitle='% Change (Referential minus Bipolar)';
 pcolorjk(binz(1:size(mb_m,2))-binsz,frx,mDiff); shading flat; ylabel('Frequency (Hz)'); xlabel('Distance (mm)'); set(gca,'fontsize',sizeoffont); colorbar; 
 title([difftitle ', ' txtyp],'fontweight','normal')
 set(gca,'ydir','normal','yscale','log','ytick',ft,'yticklabel',ftl,'xlim',xl,'ylim',fl); xlim(xldist); yline(10,'k-',.75); 
    caxdiff=clim; caxdiff=max([ abs(clim)])*[-1 1]; 
    if caxdiff(2)<100; caxdiff=[-1 1]*100; end
    clim(caxdiff); colormap(gca,cmocean('balance',40))
    %hold on; contour(binz(1:size(mb_m,2))-binsz,frx,mDiff,'k','LineWidth',0.5);
    
subplot(3,2,4)
[clusters, p_values, t_sums, permutation_distribution ] = permutest(mb_(:,:,:),mARb_(:,:,:),0); 
msig=false(nfrx,length(binz)); 
for i=1:length(clusters); 
    if p_values(i)<psig; 
        msig(clusters{i})=true; 
    end; 
end; 
mDiffcopy= mDiff;
mDiffcopy(~msig)=nan;
pcolorjk(binz(1:size(mb_m,2))-binsz,frx,mDiffcopy); shading flat; ylabel('Frequency (Hz)'); xlabel('Distance (mm)'); set(gca,'fontsize',sizeoffont); colorbar; 
 text(max(xlim)+diff(xlim)/4,mean(ylim),'% diff','fontsize',12,'rotation',90,'horizontalalignment','center')
 title([difftitle ', ' txtyp ', p<' num2str(psig)],'fontweight','normal')
 set(gca,'ydir','normal','yscale','log','ytick',ft,'yticklabel',ftl,'xlim',xl,'ylim',fl); xlim(xldist); yline(10,'k-',.75); 
    clim(caxdiff); colormap(gca,cmocean('balance',40))
    %hold on; contour(binz(1:size(mb_m,2))-binsz,frx,mDiff,'k','LineWidth',0.5);

% Histogram of number of pairs per bin
 subplot(11,2,21); histogram(make1d(Mbp_distance),[0.001 binsz:binsz:85],'facecolor',.5*[1 1 1]); set(gca,'fontsize',12); xlabel('Binned bipolar distance (mm)','fontsize',sizeoffont); 
 ylabel('Counts/bin','fontweight','normal'); axis tight; grid on; cb=colorbar; set(cb,'visible','off'); xlim(xldist); yline(10,'k-',.75); set(gca,'fontsize',sizeoffont); 
 title(pt)
 cd('~/Desktop/')
 
 savepath = '~/Desktop/bipolar_results/';
 if ~exist(savepath); mkdir(savepath); end
 saveas(gcf,[savepath pt '.png'])

%%
figure('color','w','position',[187 49 381 272]); hold on
c1=35; %using electrodes 35 and 52 for EC175 example
c2=52;
RA =squeeze(M      (c1,c1,:,:)); %channel A in referential
RB =squeeze(M      (c2,c2,:,:)); %channel B in referential
RAB=squeeze(M_averef(c1,c2,:,:)); %channe`ls A and B in referential averaged together
Rbp=squeeze(M      (c1,c2,:,:)); %channels A and B in bipolar
% transform step
RA =sqrt(RA); 
RB =sqrt(RB); 
RAB=sqrt(RAB);
Rbp=sqrt(Rbp); 
% mean step
RAm =mean(RA,2); 
RBm =mean(RB,2); 
RABm=mean(RAB,2);
Rbpm =mean(Rbp,2); 

plot(frx,(RAm)  ,'--','linewidth',1,'color',[.5 .5 .5]    ); 
plot(frx,(RBm)  ,'-.','linewidth',1,'color',[.5 .5 .5])
  %ribbons(frx,(RAB)',[0 1 1]*.75); hold on
plot(frx,(RABm),'-' ,'linewidth',2,'color',[0 0 0])
  %ribbons(frx,(Rbp)',[0 0 0]*.75); hold on
plot(frx,(Rbpm) ,'-' ,'linewidth',2,'color',[0 0 1]    )
set(gca,'xscale','log','xtick',ft,'xticklabel',ftl,'xlim',fl); grid on
legend({[num2str(35) ' Ref'],[num2str(52) ' Ref'],'Ref Ave','Bipolar'})
title(['Electrode ' num2str(c1) ' vs. ' num2str(c2) ' (' num2str(round(Mbp_distance(c1,c2),1)) ' mm)'],'fontweight','normal')
 %cd('~/Desktop/')
 saveas(gcf,[savepath pt '__ex.png']);

