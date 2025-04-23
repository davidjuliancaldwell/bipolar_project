% Script to plot patient brains and electrodes (Figure 1)
% Additionally plot stacked histograms for number of bipolar pairs
% stratified by electrode type

%%

function bipolar_electrodemaps_2025(subdural1depths2)
if ~exist('subdural1depths2','var'); subdural1depths2=1; end % 1: show lateral subdural electrdoe view, 2: show inferior view of depths only
cd /Volumes/KLEEN_DRIVE/bipolar_expedition/
fig = figure(100+subdural1depths2);
set(fig, 'Position', [100, 100, 1400, 900]);
set(gcf, 'Color', 'white');
pts = {'EC133', 'EC175', 'EC181', 'EC183', 'EC186', 'EC187', 'EC196', 'EC219', 'EC220', 'EC221', 'EC222',};%  'EC131','EC143','EC157','EC162','EC168'};



pts_names = {'Pt. 1', 'Pt. 2', 'Pt. 3', 'Pt. 4', 'Pt. 5', 'Pt. 6', 'Pt. 7', 'Pt. 8', 'Pt. 9', 'Pt. 10', 'Pt. 11'};%   'Pt. 12', 'Pt. 13', 'Pt. 14', 'Pt. 15', 'Pt. 16'};

%Brain plot positions

pos = {[0.01 0.765 0.18 0.18],
    [0.2 0.765 0.18 0.18],
    [0.39 0.765 0.18 0.18],
    [0.01 0.535 0.18 0.18],
    [0.2 0.535 0.18 0.18],
    [0.39 0.535 0.18 0.18],
    [0.01 0.305 0.18 0.18],
    [0.2 0.305 0.18 0.18],
    [0.39 0.305 0.18 0.18]
    [0.06 0.075 0.18 0.18]
    [0.25 0.075 0.18 0.18]};

%Electrode identifiers for 9 pts with depths, pulled from AN_info

depths = {[281:320],
          [289:340],
          [],
          [257:286 289:298],
          [257:266 289:318],
          [299:318],
          [273:282 289:308],
          [399:408 429:468],
          [],
          [353:382 385:404],
          [257:298]};

%Grids/strips in blue, depths in red
for i = 1:length(pos)
    subplot('Position', pos{i});
    elecsbrain(pts{i}, 0, [], [0 0 1], 'b', 0, 2.75, 2);
    elecsbrain(pts{i}, 0, depths{i}, [1 0 0], 'b', 0, 3.3, 2); 
    lightsout; 
    if subdural1depths2==1
       litebrain('l', 1);
       alpha 0.5;
       title(pts_names{i});
    elseif subdural1depths2==2
       litebrain('i', 1);
       alpha 0.3;
    end
end

plot_pos = {[0.62 0.715 0.37 0.25],
    [0.62 0.05 0.37 0.25],
    [0.62 0.38 0.37 0.25]};

for j = [1, 3, 2]
    plot_dists(j, plot_pos{j});
end

end

%%

function plot_dists(elec_type, plotposition)

    g1s2d3 = elec_type;
    pts_all = {'EC133', 'EC175', 'EC181', 'EC183', 'EC186', 'EC187', 'EC196', 'EC219', 'EC220', 'EC221', 'EC222'};
    maxs = [];

    % Define colors for each patient
    color_keys = {'EC133', 'EC175', 'EC181', 'EC183', 'EC186', 'EC187', 'EC196', 'EC219', 'EC220', 'EC221', 'EC222'};
    pts_names = {'Pt. 1', 'Pt. 2', 'Pt. 3', 'Pt. 4', 'Pt. 5', 'Pt. 6', 'Pt. 7', 'Pt. 8', 'Pt. 9', 'Pt. 10', 'Pt. 11'};

    color_values = [
        0.0, 0.4470, 0.7410;   
        0.8500, 0.3250, 0.0980; 
        0.9290, 0.6940, 0.1250; 
        0.4940, 0.1840, 0.5560;
        0.4660, 0.6740, 0.1880;
        0.3010, 0.7450, 0.9330;
        0.6350, 0.0780, 0.1840; 
        0.75, 0.75, 0;         
        0.25, 0.25, 0.25;      
        1, 0.6, 0.2;           
        0, 0.5, 0];

    color_map = containers.Map(color_keys, mat2cell(color_values, ones(1, length(color_keys)), 3));

    if g1s2d3 == 2
        load('strips_distances.mat'); % From get_elecs, all distance
        sel = [1, 2, 4, 5, 6, 7, 8, 10, 11];
        pts = pts_all(sel);
        disp('Creating strips stack...');
        elec = 'Strip';
    elseif g1s2d3 == 3
        load('depths_distances.mat');
        sel = [1:11];
        pts = pts_all(sel);
        disp('Creating depths stack...');
        elec = 'Depth';
    else
        load('grids_distances.mat');
        sel = [1, 2, 4, 5, 6, 7, 8, 10, 11];
        pts = pts_all(sel);
        disp('Creating grids stack...');
        elec = 'Grid';
    end

    subplot('Position', plotposition);
    hold on;
    grid on;

    bin_edges = 0:1:80; % Excluding distances above 80 mm
    hist_counts = zeros(length(pts), length(bin_edges) - 1);

    for i = 1:length(pts)
        pt = pts{i};
        distances = distance{sel(i)}; 
        disp(size(distances));
        
        if ~isempty(distances)
            hist_counts(i, :) = histcounts(distances, bin_edges);
        end
    end


    stacked_counts = zeros(size(hist_counts));
    h = bar(bin_edges(1:end-1), hist_counts', 'stacked', 'EdgeColor', 'none');

    for p = 1:length(pts)
        set(h(p), 'FaceColor', color_map(pts{p}));
    end

    ylabel('Counts', 'FontSize', 9);
    xlabel('Distance (mm)', 'FontSize',9);
    xlim([0 80]);
    set(gca, 'GridAlpha', 0.5);
    title([elec, ' Electrodes'], FontWeight='normal');
    
    if g1s2d3==3
        legend(pts_names, 'Location', 'Best');
        set(legend, 'String', pts_names); %de-identify
    end
    
    hold off;
    

end


