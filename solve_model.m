function [V,Policy,StationaryDist,stat_age,vf_time,time_distrib,age_time] = ...
    solve_model(par,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,...
    DiscountFactorParamNames,vfoptions,jequaloneDist,AgeWeightsParamNames,...
    simoptions,FnsToEvaluate)

disp('USING VFI TOOLKIT')

disp('Start ValueFnIter')
tic;
[V, Policy]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);
vf_time=toc;

disp('Start distribution...')
tic
StationaryDist=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightsParamNames,Policy,n_d,n_a,n_z,N_j,pi_z,Params,simoptions);
time_distrib=toc;

disp('Start AgeConditionalStats...')

tic
stat_age=LifeCycleProfiles_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions);
%[stat_age] = LifeCycleProfiles_ale_v2(pi_z,StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions);
age_time=toc;

end %end function