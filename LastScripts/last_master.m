
function last_master

% finishing analyses, outputting data for final plots

do_fig2_grids = true;
do_fig4 = false; % set to true if havent executed yet
do_fig3 = true;
do_supp1 = true;

if do_fig2_grids

    bipolarexpedition_Linear_2025Grid(1,2); % save grid linear (sqrt)

    % Note --> can recompute supp fig 1 by simply accessing dk_lin variable
end

if do_fig4
    loopbipolarexpedition(2); % Strips Fig 4
    loopbipolarexpedition(3); % Depths Fig 4
end

% Fig 3

if do_fig3
    pts = {'EC175', 'EC183'};
    for i = 1:size(pts,2)
        [pt, binz, toplot, frx, binsz, Mbp_distance, cm_distance] = fig3_EachVsAll(pts{i});
        save(['/Users/jonathankleen/Desktop/bipolar_results/' pts{i} '_fig3_data.mat'], 'pt', 'binz', 'toplot', 'frx', 'binsz', 'Mbp_distance', 'cm_distance');
        disp(['Finished: ' pts{i}]);
    end
end

% Supp 1 --> Recalculate for freq bands

% Edit DK: can recompute just from dk_lin (mean(mean(sqrt(power)))

end
