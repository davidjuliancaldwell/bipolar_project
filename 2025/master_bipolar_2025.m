%% Master script for running bipolar analysis
% Includes figures

function master_bipolar_2025

% Figure 1, electrode maps
bipolar_electrodemaps_2025;

% Figure 2, linear analysis
bipolarexpedition_Linear_2025 % View documentation with saving/plotting


% Figure 3, Each Vs All Analysis (Pt. 2, 4)
fig3_call; % View documentation with saving/plotting

% Fig 4, Relative Change from Referential
pt = 'EC175'; % or 183
EachVsAll_cleaned2025(pt);


% Figure 5, STG, task-based (Pt. 2, 4)
fig5_out;

% Figure 6
high_density_ecog_script;
postHocStats;

% Supplementary Figure 2
byconds;


end