%% Master script for running bipolar analysis
% Includes figures

function master_bipolar_2025

none1sqrt2log3=2; % 1: no transform, 2: square root, 3: log

% %Figure 1, electrode maps
% bipolar_electrodemaps_2025(1);
% bipolar_electrodemaps_2025(2);



%Figure 2, linear analysis

% View documentation with saving/plotting (1st argument is g1s2d3)
% bipolarexpedition_Linear_2025(1,none1sqrt2log3) % grids
% bipolarexpedition_Linear_2025(2,none1sqrt2log3) % strips
% bipolarexpedition_Linear_2025(3,none1sqrt2log3) % depths



%Figure 3, Each Vs All Analysis (Pt. 2, 4)
fig3_call(none1sqrt2log3); % View documentation with saving/plotting



%Fig 4, Relative Change from Referential
pt = 'EC175'; % or 183
EachVsAll_cleaned2025(pt,none1sqrt2log3);



%Figure 5, STG, task-based (Pt. 2, 4)
fig5_out;



%Figure 6
high_density_ecog_script;
postHocStats;


end