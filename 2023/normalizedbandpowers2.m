%% NORMALIZING BY FREQ BAND AND BIPOLAR DISTANCE
%For High Density Grids Only

bandpowers = load('/scratch/bipolar_expedition/normalizedbandpowers.mat');
pts = {'EC133','EC175','EC183','EC186','EC187','EC196','EC219','EC221','EC222'};

delta = bandpowers.delta;
theta = bandpowers.theta;
alpha = bandpowers.alpha;
beta = bandpowers.beta;
gamma = bandpowers.gamma;
highgamma = bandpowers.highgamma;

colors = {[0, 1.0, 0], [0, 0.5, 1.0], [1.0, 0, 0], [0.75, 0, 0.75], [1.0, 0.65, 0], [0, 1.0, 1.0], [0, 0, 0.5], [1.0, 0, 1.0], [0, 0.75, 0.25]};

aggregate = {delta, theta, alpha, beta, gamma, highgamma};
subtitles = {'Delta (2-4Hz)', 'Theta (4-8Hz)', 'Alpha (8-13Hz)', 'Beta (13-25Hz)', 'Gamma (25-50Hz)', 'High Gamma (50-1200Hz)'};
ticklabels = {'Referential', '4 mm', '8 mm', '12 mm', '16 mm', '20 mm'};

fig = figure();
set(gcf, 'Position', [100, 100, 1400, 750]); 

handles = gobjects(1, length(pts));

for i = 1:length(subtitles)
    subplot(3, 3, i) 
    for j = 1:length(pts)
         h = plot(aggregate{i}(j,:), 'Color', colors{j}, 'LineWidth', 2, 'Marker', 'o', 'MarkerSize', 4, 'MarkerFaceColor', colors{j});
         handles(j) = h; 
         hold on;
    end
    title(subtitles{i}, 'FontSize', 14);
    ylim([-1.5 1.0]);
    yticks(-1.5:0.5:1.0);
    xticks(1:6);
    xticklabels(ticklabels);
    xlabel('Bipolar distance');
    ylabel('Mean log power');
    grid on
end


subplot(3, 3, [7 9]); 
axis off; 
legend(handles, pts, 'Location', 'north', 'Orientation', 'vertical', 'EdgeColor', 'none', 'FontSize', 12);



%%

% RUN RMANOVA on desired band

function rm_run(band);

data = band;
data = data(:, 2:end);
T = array2table(data, 'VariableNames', {'_4_mm', '_8_mm', '_12_mm', '_16_mm', '_20_mm'});
Meas = table([1 2 3 4 5]', 'VariableNames', {'Groups'});
rm = fitrm(T, '_4_mm-_20_mm ~ 1', 'WithinDesign', Meas);
ranovatbl = ranova(rm);
disp(ranovatbl);
return

end

























