clear
clc
close all

% Run main with do_toolkit = 1, 
% save results_toolkit.mat
% Then run main with do_toolkit = 0
% save results_ale.mat
% Finally, run this script to compare results


load results_toolkit.mat

% Convert to (a,h,z,j)

V_t = V;
Policy_t = reshape(Policy,[2,n_a(1),n_a(2),prod(n_z),N_j]);
StationaryDist_t = reshape(StationaryDist,[n_a(1),n_a(2),prod(n_z),N_j]);
l_coh_toolkit = stat_age.l_coh.Mean';

load results_ale.mat

V_ale = V;
Policy_ale = Policy;
StationaryDist_ale = StationaryDist;
l_coh_ale = stat_age.l_coh;

max(abs(Policy_t-Policy_ale),[],"all")
max(abs(StationaryDist_t-StationaryDist_ale),[],"all")

max(abs(l_coh_toolkit-l_coh_ale))

figure
plot(1:Params.Jr-1,l_coh_toolkit(1:Params.Jr-1))
hold on
plot(1:Params.Jr-1,l_coh_ale(1:Params.Jr-1))
legend('toolkit','my codes')

