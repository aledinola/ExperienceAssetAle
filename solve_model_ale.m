function [V,Policy,StationaryDist,stat_age,vf_time,time_distrib,age_time] = ...
    solve_model_ale(par,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,...
    DiscountFactorParamNames,vfoptions,jequaloneDist,AgeWeightsParamNames,...
    simoptions,FnsToEvaluate)

z_grid = gather(z_grid);
pi_z   = gather(pi_z);

disp('Start ValueFnIter')
tic;
% See _v2 for the vectorized version
if par.do_vfi_vec==0
    % on CPU
    [V,Policy] = ValueFnIter_Case1_FHorz_ale(n_d,n_a,n_z,N_j, d_grid, a_grid, z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames,[], vfoptions,par);
elseif par.do_vfi_vec==1
    % GPU, low memory usage 
    [V,Policy] = ValueFnIter_Case1_FHorz_ale_v2(n_d,n_a,n_z,N_j, d_grid, a_grid, z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames,[], vfoptions,par);
elseif par.do_vfi_vec==2
    % GPU, more vectorized, high memory usage
    [V,Policy] = ValueFnIter_Case1_FHorz_ale_v3(n_d,n_a,n_z,N_j, d_grid, a_grid, z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames,[], vfoptions,par);
else
    error('do_vfi_vec out of bounds')
end
vf_time=toc;

% Distribution and life-cycle profiles do not use GPU
V = gather(V);
Policy = gather(Policy);

disp('Start distribution...')
tic
StationaryDist=StationaryDist_FHorz_Case1_ale(par,z_grid,jequaloneDist,AgeWeightsParamNames,Policy,n_d,n_a,n_z,N_j,pi_z,Params,simoptions);
time_distrib=toc;

disp('Start AgeConditionalStats...')

tic;
[stat_age] = LifeCycleProfiles_ale(pi_z,StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions);
age_time=toc;

end %end function