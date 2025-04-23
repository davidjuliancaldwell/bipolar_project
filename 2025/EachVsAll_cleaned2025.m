% BIPOLAR PAIR ANALYSIS: EACH VS. ALL to make figures
function [mDiff, mb_m, mARb_m,Mbp_distance] = EachVsAll_cleaned2025(pt,g1s2d3,none1sqrt2log3);
%function Mbp_distance = devon_EachVsAll_cleaned(pt)

%if ~exist('pt','var')||isempty(pt); pt='EC175'; end %pt='EC175'; % EC175 and EC183 both have intact 16x16 square grids (channel #s 1:256)
%if ~exist('nchtocheck','var')||isempty(nchtocheck); nchtocheck=128*2; end
%if ~exist('windowstocheck','var')||isempty(windowstocheck); windowstocheck=250; end %each window is 1 second of data (non-overlapping)
% pt='EC175'; %pt='EC175'; % EC175 and EC183 both have intact 16x16 square grids (channel #s 1:256)
windowstocheck=250; %each window is 1 second of data (non-overlapping)
% g1s2d3=1; % use either grids (1) or strips (2) or depths (3) but not the others
% none1sqrt2log3 values for the pre-calculation transform: 1: no transform, 2: square root, 3: natural log
    if ~exist('none1sqrt2log3','var'); none1sqrt2log3=2; msgbox(['FYI, defaulting to ' Txform{none1sqrt2log3} ', use as input next time']); end
doanglerange=0;
recordings = [];

cm=cool(6); cm(1,:)=[0 0 0]; 

data_root = getenv("KLEEN_DATA");
datadir = fullfile(data_root, 'bipolar_expedition', 'baseline-high-density-data');
tag_spikes_path = fullfile(data_root, 'bipolar_expedition', 'taggedspikes_April2022.mat');
load(tag_spikes_path);

sfx=512;
frxrange=[2 200]; %frequency range to examine
ft=[2 5 10 20 50 100 200]; ftl=cellstr(num2str(ft')); %frequency labels for plots
 
u=dir(datadir); uptbl={}; for i=1:length(u); uname=u(i).name; uptbl{i,1}=uname(1:end-28); end; uptbl(1:2)=[]; clear i u uname

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
    ptpath = fullfile(datadir, [uptbl{ptbl(b)} '_baselineWindows_fromraw.mat']);
    load(ptpath);
    recordings = [recordings size(nonspks_windows, 2)];
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

an_electrode_info_path = fullfile(data_root, 'bipolar_expedition', 'AN_ElectrodeInfoTDT.xlsx');
[bpN,bpT]=xlsread(an_electrode_info_path, pts{p});
[em,eleclabels,anatomy]=getelecs(pts{p},2);

cm=cool(6); cm(1,:)=[0 0 0]; 

datadir = fullfile(datadir, 'bandpassfiltered');

if g1s2d3 == 1
    use_ch = find(strcmpi('grid', anatomy(:,3)) | strcmpi('minigrid', anatomy(:,3)));
    binsz=2; % bin size in mm
    xldist = [0 60];
elseif g1s2d3 == 2
    use_ch = find(strcmpi('strip', anatomy(:,3)));
    binsz = 4;
    xldist = [0 60];
else
    use_ch = find(strcmpi('depth', anatomy(:,3)));
    binsz = 4; %binsz=2;
    xldist = [0 60];
end

badchI=isnan(mean(mean(d,1),3))'; %any channel with nans in any window will be a "bad channel"
badchI(badchidx{p})=1; %follow up with all marked as bad channels in original preprocessing
x=[]; xbch=true(size(badchI)); for i=1:size(bpT,1); x=[x bpN(i,1):bpN(i,2)]; end; xbch(x)=false; 
badchI(xbch)=1; %find and nan empty channels unaccounted for by component rows
d(:,badchI,:)=nan;
okc=~badchI; clear x xbch

windowstocheck=min([windowstocheck size(d,3)]);
windowstocheck=1:windowstocheck; %convert to a vector of windows, 1:X

% ALL PAIRS (each vs. all others) analysis and example plot
 d=d(:,:,windowstocheck); clear Straces_allch; %free up RAM by getting rid of whatever won't be used (only using first ___ number of windows)
                 % ***opportunity here to select speech or stim windows
 d = d(:, use_ch, :);
 nchtocheck = size(d, 2);

 
 
 % NEW DISTANCE CALCULATIONS
 %{
addons = [];
for i = 1:nchtocheck
    val = use_ch(i)-i;
    addons = [addons val];
end
 

% euclidean distance and angle (in sagittal plane) for each pair

for c1=1:nchtocheck; 
     for c2=1:nchtocheck; 
         Mbp_distance(c1,c2)=distance3D(em(c1+addons(c1),:),em(c2+addons(c1),:)); 
         %get angle in sagittal plane for each pair (all L-side pts so no need to flip anyone)
          % this is from the inverse tangent of the vertical (superior/inferior) axis difference divided by
          % horizontal (anterior/posterior) axis difference of the two electrodes
          if c1~=c2; Mbp_angle   (c1,c2)=atan2((em(c2,3)-em(c1,3)),(em(c2,2)-em(c1,2))); end
     end; 
 end
    rmv=Mbp_distance==0; %remove spurious zero distance values
    Mbp_distance(rmv)=nan; 
    Mbp_angle   (rmv)=nan; 
    %referential signal denoted as "0" distance for coding/indexing purposes below
    for c1=1:nchtocheck; 
        Mbp_distance(c1,c1)=0; 
        Mbp_angle   (c1,c1)=nan; % angle will still be nans (to avoid confusion with 0 degrees)
    end 
 % remove mirror image in confusion matrix
 for c2=1:nchtocheck; 
     for c1=c2+1:nchtocheck; 
         M(c1,c2,:)=nan; 
         Mbp_distance(c1,c2)=nan; 
         Mbp_angle   (c1,c2)=nan;
     end; 
 end
 Mbp_distance(rmv)=nan; 
 %}

%%
 % Here

 [M,Mrefave,Mbp_distance,frx,~,Mbp_angle]=bpspectra_EachVsAll_2023(d,sfx,frxrange,em,nchtocheck);

% [M,Mrefave,frx,~]=dmod_bpspectra_EachVsAll_2023(d,sfx,frxrange,em,nchtocheck);

     % M is bipolarchannel1 X bipolarchannel2 X frx X 1secwindow
 nfrx=length(frx);
  

% ANGLE RANGE: if desired, subselect bipolar angle within a given range

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

% Unstack and line up into 2D matrix to create channels^2 X frequencies for easier indexing-->binning
mb=[]; 
mARb=[]; 
binz=[-1 0:binsz:85]; % can change binsz PRN at this point to test different bin resolutions on the plots below
clear mx my mz  
nbinz=length(binz)-1;
 parfor w=windowstocheck; disp(num2str(w)); % parfor here to run the windows
    Mflat=[];
      Mflatbp_distance=[];
      %Mflatbp_angle=[];
    Mrefaveflat=[];
    for c2=1:nchtocheck; 
        Mflat=[Mflat; squeeze(M(:,c2,:,w))];               
        Mflatbp_distance=[Mflatbp_distance; Mbp_distance(:,c2)]; % corresponding distance index
        %Mflatbp_angle   =[Mflatbp_angle;    Mbp_angle(:,c2)]; % corresponding angle index
        Mrefaveflat=[Mrefaveflat; squeeze(Mrefave(:,c2,:,w))];              
    end
    
    %bin by distance and take the mean, creating: frequency X binned distance
       % first for referential (distance = 0)
    for i=1:nbinz % binz will be including >lower bound and up to and including (<=) upper bound
        mb  (:,i,w)=nanmean(Mflat      (Mflatbp_distance>binz(i) & Mflatbp_distance<=binz(i+1),:),1); 
        mARb(:,i,w)=nanmean(Mrefaveflat(Mflatbp_distance>binz(i) & Mflatbp_distance<=binz(i+1),:),1);
    end
    
    disp([num2str(round(windowstocheck(w)/windowstocheck(end)*100,1)) '% of windows'])
end; disp('Done')

binz(1)=[]; 

%% zscore the log transformed power according to frequency

mb__m=squeeze(mean(sqrt(mb),3)); 
mARb__m=squeeze(mean(sqrt(mARb),3)); 
 nfrx=length(frx);
 mb__m_z=mb__m; mARb__m_z=mARb__m;
 for i=1:nfrx; nns=~isnan(mb__m(i,:)); 
     mb__m_z(i,nns)=zscore((mb__m(i,nns))); 
     mARb__m_z(i,nns)=zscore((mARb__m(i,nns))); 
 end; 


 %% plot log-transformed version
 
 
 figure('color','w','position',[[54 223 1907 1102]]); colormap(parula); 
 toplot=mean((mb__m_z),3);  % mb or mb_z
 cm_distance=flipud(cmocean('deep',1+ceil(max(max(Mbp_distance)))));

% Plot confusion matrix of bipolar distances between all pairs
 subplot(8,3,1:3:7); pcolorjk(Mbp_distance); shf; hold on; %plot([1 1; 1 256]',[1 256;256 256]','k-','linewidth',1); plot([1 1 1 128 256],[1 128 256 256 256],'k.'); %plot([1 128 256 256 256]/2,[1 1 1 128 256]/2,'k.'); 
    axis equal off; set(gca,'ydir','normal','xdir','reverse'); colormap(gca,cm_distance); colorbar('fontsize',12); title('Bipolar pair distance (mm)','fontsize',14,'fontweight','normal')
    clim([0 80]); 
% Plot confusion matrix of bipolar distances between all pairs
 subplot(8,3,2:3:8); pcolorjk(rad2deg(Mbp_angle)); shf; hold on; %plot([1 1; 1 256]',[1 256;256 256]','k-','linewidth',1); plot([1 1 1 128 256],[1 128 256 256 256],'k.'); %plot([1 128 256 256 256]/2,[1 1 1 128 256]/2,'k.'); 
    axis equal off; set(gca,'ydir','normal','xdir','reverse'); colormap(gca,hsv(360)); clim([-180 180]); colorbar('fontsize',12); title('Bipolar pair angle (deg)','fontsize',14,'fontweight','normal')

% Histogram of number of pairs per bin
 subplot(8,2,7); histogram(make1d(Mbp_distance),0:binsz:85,'facecolor',.5*[1 1 1]); set(gca,'fontsize',12); xlabel('Binned bipolar distance (mm)','fontsize',14); 
 ylabel('Counts/bin','fontweight','normal'); axis tight; grid on; cb=colorbar; set(cb,'visible','off'); xlim(xldist)

%saveas(gcf, '/home/devkrish/bipolar_project/2023/output/LogTransEachVAll.png', 'png');
 

% Histogram of angles

if doanglerange
    subplot(8,6,22); angs=make1d(Mbp_angle); nns=~isnan(angs);
 polarhistogram(angs(nns),-pi:2*pi*(5/360):pi,'facecolor',.5*[1 1 1]); title('Bipolar angle distribution (5^o steps)','fontsize',14,'fontweight','normal'); 
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
 pcolorjk(binz(1:size(toplot,2)),frx,toplot); shading flat; set(gca,'ydir','normal'); ylabel('Frequency (Hz)'); xlabel('Distance (mm)'); set(gca,'fontsize',14); colorbar; 
 title({'ln(power), z-scored by frequency',''},'fontweight','normal')
 set(gca,'yscale','log','ytick',ft,'yticklabel',ftl); xlim(xldist); clim([-1 1]*(max(abs(clim)))); colormap(gca,cmocean('balance')); %cm=cmocean('balance',100); colormap(gca,cm(:,[2 1 3]))
 clim([-4.5 4.5]);


 subplot(2,2,4)
 [mx,my]=meshgrid(binz(1:size(toplot,2))+binsz/2,frx);
 surf(mx,my,(toplot)); xlabel('Bipolar distance (mm)','fontsize',14); ylabel('Frequency (Hz)','fontsize',14); zlabel({'ln(power),','z-scored by frequency'},'fontsize',14)
 view(140,25); set(gca,'xdir','reverse','ydir','reverse','ytick',ft,'yticklabel',ftl,'fontsize',14); 
 colormap(gca,cmocean('balance')); clim([-1 1]*(max(abs(clim)))); 
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
% will also add a line at 10mm for visualization of this common clinical inter-electrode distance

% transform power before calculating perecent change
% but first, implement a transform, 
%       keep in mind averaging and discerning % change are infleunced 
%       strongly by which transform is used
% none1sqrt2log3 values: 1: no transform, 2: square root, 3: natural log
mb_=mb; mARb_=mARb; % make a copy... then transform the copy
if none1sqrt2log3==1;     txtyp='raw'; % no transform, power in raw form
    mb_(~isnan(mb_))      =    (mb_(~isnan(mb_)));
    mARb_(~isnan(mARb_))  =    (mARb_(~isnan(mARb_)));
elseif none1sqrt2log3==2; txtyp='square root'; % square root transform
    mb_(~isnan(mb_))      =sqrt(mb_(~isnan(mb_)));
    mARb_(~isnan(mARb_))  =sqrt(mARb_(~isnan(mARb_)));
elseif none1sqrt2log3==3; txtyp='natural log'; % log transform --> problematic because taking % change of certain small negative values creates extreme values
    mb_(~isnan(mb_))      =log (mb_(~isnan(mb_))+1);
    mARb_(~isnan(mARb_))  =log (mARb_(~isnan(mARb_))+1);
end
    
% mean across windows
mb_m=squeeze(mean(mb_,3)); 
mARb_m=squeeze(mean(mARb_,3)); 
mDiff=((mb_m-mARb_m)./mARb_m)*100;


psig=0.05;

xl=[binsz*2 binz(end)]; fl=frx([1 end]);

figure('color','w','position',[315 1 1486 1046]); colormap(jet); clear cax; %
 subplot(3,2,1)
 pcolorjk(binz(1:size(mb_m,2))-binsz,frx,mARb_m); shading flat; ylabel('Frequency (Hz)'); xlabel('Distance (mm)'); set(gca,'fontsize',14); colorbar; 
 text(max(xlim)+diff(xlim)/4,mean(ylim),'(power)','fontsize',12,'rotation',90,'horizontalalignment','center')
 title({'Referential',['(average of pair)' ', ' txtyp]},'fontweight','normal')
 set(gca,'ydir','normal','yscale','log','ytick',ft,'yticklabel',ftl,'xlim',xl,'ylim',fl);
 cax(1,:)=clim;

 subplot(3,2,2)
 pcolorjk(binz(1:size(mb_m,2))-binsz,frx,mb_m); shading flat; ylabel('Frequency (Hz)'); xlabel('Distance (mm)'); set(gca,'fontsize',14); colorbar; 
 text(max(xlim)+diff(xlim)/4,mean(ylim),'(power)','fontsize',12,'rotation',90,'horizontalalignment','center')
 title(['Bipolar' ', ' txtyp],'fontweight','normal')
 set(gca,'ydir','normal','yscale','log','ytick',ft,'yticklabel',ftl,'xlim',xl,'ylim',fl);
 cax(2,:)=clim;

 cax=(([min(cax(:,1),[],1) max(cax(:,2),[],1)])); % get extreme min and max of the two conditions
 subplot(3,2,1);   clim(cax); xlim(xldist); yline(10,'k-',.75); 
 subplot(3,2,2);   clim(cax); xlim(xldist); yline(10,'k-',.75); 

 subplot(3,2,3)
 %mDiff=((mb_m-mARb_m)./mARb_m)*100; 
 difftitle='% Change (Referential minus Bipolar)';
 pcolorjk(binz(1:size(mb_m,2))-binsz,frx,mDiff); shading flat; ylabel('Frequency (Hz)'); xlabel('Distance (mm)'); set(gca,'fontsize',14); colorbar; 
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
pcolorjk(binz(1:size(mb_m,2))-binsz,frx,mDiffcopy); shading flat; ylabel('Frequency (Hz)'); xlabel('Distance (mm)'); set(gca,'fontsize',14); colorbar; 
 text(max(xlim)+diff(xlim)/4,mean(ylim),'% diff','fontsize',12,'rotation',90,'horizontalalignment','center')
 title([difftitle ', ' txtyp ', p<' num2str(psig)],'fontweight','normal')
 set(gca,'ydir','normal','yscale','log','ytick',ft,'yticklabel',ftl,'xlim',xl,'ylim',fl); xlim(xldist); yline(10,'k-',.75); 
    clim(caxdiff); colormap(gca,cmocean('balance',40))
    %hold on; contour(binz(1:size(mb_m,2))-binsz,frx,mDiff,'k','LineWidth',0.5);

% Histogram of number of pairs per bin
 subplot(11,2,21); histogram(make1d(Mbp_distance),[0.001 binsz:binsz:85],'facecolor',.5*[1 1 1]); set(gca,'fontsize',12); xlabel('Binned bipolar distance (mm)','fontsize',14); 
 ylabel('Counts/bin','fontweight','normal'); axis tight; grid on; cb=colorbar; set(cb,'visible','off'); xlim(xldist); yline(10,'k-',.75); set(gca,'fontsize',14); 
 title(pt)

 %cd('~/Desktop/')
 %saveas(gcf,['/home/devkrish/bipolar_project/2023/output/' pt '.png'])
%}
 
%% 


figure('color','w','position',[104 374 381 431]); hold on

if g1s2d3 == 1 %from EC175
    rename=0;
    c1=35;
    c2=52;
elseif g1s2d3 == 3 %from EC181
    rename=1;
    c1=find(use_ch==65);
    c2=find(use_ch==66);
else
    rename=2;
    c1=11;
    c2=12;
end


RA =squeeze(M      (c1,c1,:,:)); %channel A in referential
RB =squeeze(M      (c2,c2,:,:)); %channel B in referential
RAB=squeeze(Mrefave(c1,c2,:,:)); %channe`ls A and B in referential averaged together
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


plot(frx,(RAm)  ,'--','linewidth',3,'color',[0 0 1]    ); 
plot(frx,(RBm)  ,'--','linewidth',3,'color',[0 1 0]*.75)
  %ribbons(frx,(RAB)',[0 1 1]*.75); hold on
plot(frx,(RABm),'-' ,'linewidth',2,'color',[0 1 1]*.75)
  %ribbons(frx,(Rbp)',[0 0 0]*.75); hold on
plot(frx,(Rbpm) ,'-' ,'linewidth',2,'color',[0 0 0]    )
set(gca,'xscale','log','xtick',ft,'xticklabel',ftl,'xlim',fl); grid on


if rename==1
    title(['Ch ' num2str(65) ' to ' num2str(66) ' -- ' num2str(Mbp_distance(c1,c2)) 'mm'])
elseif rename==0
    title(['Ch ' num2str(c1) ' to ' num2str(c2) ' -- ' num2str(Mbp_distance(c1,c2)) 'mm'])
else
    title(['Ch ' num2str(309) ' to ' num2str(310) ' -- ' num2str(Mbp_distance(c1,c2)) 'mm'])
end

 %cd('~/Desktop/')
 savepath = fullfile('~', 'Desktop', [pt '__ex.png']);
 %saveas(gcf,savepath);

%}

%% FIGURE OUTPUT

% Chs Plot


fig = figure(4); % figure 4
set(fig, 'Position', [100, 100, 1200, 900]);

subplot1 = subplot('Position', [0.0832, 0.7277, 0.2509, 0.2344]);
plot(frx,(RAm)  ,'-','linewidth',2.5,'color',[0.6 0.6 0.6], 'DisplayName', [num2str(c1) ' Ref']);  
hold on;

plot(frx,(RBm)  ,'-','linewidth',2.5,'color',[0.6 0.6 0.6], 'DisplayName', [num2str(c2) ' Ref']);
hold on;

plot(frx,(RABm),'-' ,'linewidth',3,'color',[0 0 0], 'DisplayName', 'Ref Ave');
hold on;

plot(frx,(Rbpm) ,'-' ,'linewidth',3,'color',[0 0 1], 'DisplayName', 'Bipolar');
set(gca,'xscale','log','xtick',ft,'xticklabel',ftl,'xlim',fl, 'GridAlpha', 0.4, 'MinorGridAlpha', 0.55); 
grid on; 

if rename==1
    title(['Ch ' num2str(65) ' to ' num2str(66) ' -- ' num2str(Mbp_distance(c1,c2)) 'mm'])
elseif rename==0
    title(['Ch ' num2str(c1) ' to ' num2str(c2) ' -- ' num2str(Mbp_distance(c1,c2)) 'mm'])
else
    title(['Ch ' num2str(299) ' to ' num2str(300) ' -- ' num2str(Mbp_distance(c1,c2)) 'mm'])
end

hLegend = legend('show');
hLegend.ItemTokenSize = [15, 18];

subplot2 = subplot('Position', [0.0649, 0.6202, 0.2834, 0.0605]); 
histogram(make1d(Mbp_distance),[0.001 binsz:binsz:85],'facecolor',.5*[1 1 1]); set(subplot2,'fontsize',12); xlabel('Binned bipolar distance (mm)','fontsize',9); 
ylabel('Counts/bin','fontweight','normal'); axis tight; grid on; xlim(xldist); xline(10,'k-'); 
set(subplot2,'fontsize',9); 
title(pt);
colormap(jet); 
clear cax;
position1 = [0.0584, 0.3356, 0.32, 0.2177]; 
position2 = [0.0584, 0.0512, 0.32, 0.2177]; 

subplot3 = subplot('Position', position1);
pcolor(binz(1:size(mb_m,2))-binsz, frx, mARb_m); shading flat; 
ylabel('Frequency (Hz)'); xlabel('Distance (mm)'); 
set(gca,'fontsize',9); colorbar; 
%text(max(xlim)+diff(xlim)/4, mean(ylim), '(power)', 'fontsize',12, 'rotation',90, 'horizontalalignment','center')
title({'Referential Average, ln(Power)'})
set(gca,'ydir','normal', 'yscale','log', 'ytick',ft, 'yticklabel',ftl, 'xlim',xl, 'ylim',fl);
cax(1,:) = caxis; % Store color axis limits

subplot4 = subplot('Position', position2);
pcolor(binz(1:size(mb_m,2))-binsz, frx, mb_m); shading flat; 
ylabel('Frequency (Hz)'); xlabel('Distance (mm)'); 
set(gca,'fontsize',9); colorbar; 
%text(max(xlim)+diff(xlim)/4, mean(ylim), '(power)', 'fontsize',12, 'rotation',90, 'horizontalalignment','center')
title(['Bipolar, ln(Power)'])
set(gca,'ydir','normal', 'yscale','log', 'ytick',ft, 'yticklabel',ftl, 'xlim',xl, 'ylim',fl);
cax(2,:) = caxis; 

combined_cax = [min(cax(:,1)), max(cax(:,2))];
axes(subplot3); caxis(combined_cax);
xlim(xldist); xline(10,'k-'); 
axes(subplot4); caxis(combined_cax);
xlim(xldist); xline(10,'k-'); 

% each pt plot

pts={'EC133','EC175','EC181','EC183','EC186','EC187','EC196','EC219','EC220','EC221','EC222'};


ft=[2 5 10 20 50 100 200]; ftl=cellstr(num2str(ft'));

if g1s2d3==1
    pos_ECs = {
        'EC133', [0.41, 0.8011, 0.165, 0.1522],
        'EC175', [0.6095, 0.8011, 0.165, 0.1522], 
        'EC181'  NaN(1,4),
        'EC183', [0.809, 0.8011, 0.165, 0.1522], 
        'EC186', [0.41, 0.6121, 0.165, 0.1522], 
        'EC187', [0.6095, 0.6121, 0.165, 0.1522], 
        'EC196', [0.809, 0.6121, 0.165, 0.1522],
        'EC219', [0.41, 0.4231, 0.165, 0.1522], 
        'EC220'  NaN(1,4),
        'EC221', [0.41, 0.2342, 0.165, 0.1522], 
        'EC222', [0.41, 0.0447, 0.165, 0.1522]};
    exclude = [3, 9];
    load('m_EachVsAll.mat');
    load('frx_binz.mat');
    limptplots = 60;

elseif g1s2d3==3 %depths
    pos_ECs = {
        'EC133', [0.41, 0.8011, 0.165, 0.1522],
        'EC175', [0.6095, 0.8011, 0.165, 0.1522], 
        'EC181'  [0.809, 0.8011, 0.165, 0.1522],
        'EC183', NaN(1,4),  
        'EC186', [0.41, 0.6121, 0.165, 0.1522], 
        'EC187', [0.6095, 0.6121, 0.165, 0.1522], 
        'EC196', NaN(1,4),
        'EC219', [0.809, 0.6121, 0.165, 0.1522],
        'EC220', [0.41, 0.4231, 0.165, 0.1522], 
        'EC221', [0.41, 0.2342, 0.165, 0.1522], 
        'EC222', [0.41, 0.0447, 0.165, 0.1522]};
    exclude = [4, 7];
    load('m_EachVsAll_depths_bins4.mat');
    load('frx_binz_depths.mat');
    limptplots=60;
else
    pos_ECs = {
        'EC133', [0.41, 0.8011, 0.165, 0.1522],
        'EC175', [0.6095, 0.8011, 0.165, 0.1522], 
        'EC181'  NaN(1,4),
        'EC183', [0.809, 0.8011, 0.165, 0.1522], 
        'EC186', [0.41, 0.6121, 0.165, 0.1522], 
        'EC187', [0.6095, 0.6121, 0.165, 0.1522], 
        'EC196', [0.809, 0.6121, 0.165, 0.1522],
        'EC219', [0.41, 0.4231, 0.165, 0.1522], 
        'EC220'  NaN(1,4),
        'EC221', [0.41, 0.2342, 0.165, 0.1522], 
        'EC222', [0.41, 0.0447, 0.165, 0.1522]};
    exclude = [3, 9];
    load('m_EachVsAll_strips_bins4.mat');
    load('frx_binz_depths.mat');
    limptplots = 60;
end

subplotHandles = zeros(1, size(pos_ECs, 1));

for i = 1:length(pts)
    if isnan(pos_ECs{i,2})
        continue;
    end
    subplotHandles(i) = axes('Position', [pos_ECs{i,2}]);
    pcolorjk(binz_re,frx_re,squeeze(mDiff_re(i,:,:)));
    cax(i,:)=clim;
    title(pos_ECs{i});
end

cm = cmocean('balance',24*2);

for j = 1:length(subplotHandles)
    if subplotHandles(j) == 0
        continue
    end
    subplot(subplotHandles(j));
    cax=[-60 60];
end


for i=1:11
    if i==exclude(1) || i==exclude(2); continue; end
    subplot(subplotHandles(i)); 
    colormap(subplotHandles(i), cm);
    clim(cax); shading flat; set(gca,'yscale','log','ytick',ft,'yticklabel',ftl); %colorbar; 
    xlim([4.00 limptplots]); 
    xticks([10 20 30 40 50 60 70])
end

% all plot

if rename==1
    load('m_EachVsAll_depths_bins2.mat');
    load('frx_binz.mat');
end

if rename==2
    load('m_EachVsAll_strips_bins2.mat');
    load('frx_binz.mat');
end

all_meanplot = subplot('Position', [0.62, 0.0678, 0.363, 0.4844]);
pcolorjk(binz_re,frx_re,squeeze(nanmean(mDiff_re,1))); shading flat; title('All (mean)'); %nanmean vs mean
colormap(all_meanplot, cm);
clim(cax); shading flat; set(gca,'yscale','log','ytick',ft,'yticklabel',ftl); colorbar; xlim([4.00 20]);
ylabel('Frequency (Hz)'); xlabel('Distance (mm)'); xticks([4 6 8 10 12 14 16 18 20]);
set(gca,'fontsize',9);

%}

%% Saving

%{
%savepath = fullfile('~', 'Desktop', 'Fig4_strips_new.png');
%saveas(gcf, savepath);

%} 
