function StationaryDist=StationaryDist_FHorz_Case1_ale(par,z_grid,jequaloneDist,AgeWeightParamNames,Policy,n_d,n_a,n_z,N_j,pi_z,Params,simoptions)

verbose = simoptions.verbose;
d_grid = simoptions.d_grid; % because we are using an experience asset
a_grid = simoptions.a_grid; % because we are using an experience asset

% NS_m = n_z(1);
% NS_f = n_z(1);
% NP   = n_z(3);

[n_a,n_h,a_grid,h_grid,z_grid,n_z] = prepare_grids(n_a,a_grid,z_grid,n_z,pi_z);

% State variables:
% (a,h,eta_m,eta_f,theta,j)-->(a,h,z,j)
%jequaloneDist(1,1,eta_m_init,eta_f_init,ip)

% Policy(1,a,h,eta_m,eta_f,theta,j) for d, decision variable
% Policy(2,a,h,eta_m,eta_f,theta,j) for a', future endog state

%% Set parameters that do not depend on age
xi_1  = Params.xi_1;
xi_2  = Params.xi_2;
del_h = Params.del_h;
h_l   = Params.h_l;
%eta_m_init = par.is_initial_m;
%eta_f_init = par.is_initial_f;

%% Set initial distribution at age 1
mu = zeros(n_a,n_h,n_z,N_j);
% Go from (a,h,eta_m,eta_f,theta)-->(a,h,z)
% Note that n_z = NS_m*NS_f*NP
jequaloneDist = reshape(jequaloneDist,[n_a,n_h,n_z]);
mu(:,:,:,1) = jequaloneDist;

%% Forward iteration

for j=2:N_j

    % Set age-dependent parameters
    age_j = Params.agej(j);

    for z_c=1:n_z
    for h_c=1:n_h
    for a_c=1:n_a
        % Next-period a'
        aprime = Policy(2,a_c,h_c,z_c,j-1);
        % Decision variable d
        d_c = Policy(1,a_c,h_c,z_c,j-1);
        % Next-period h'
        d_val = d_grid(d_c);
        h_val = h_grid(h_c);
        hprime_val = f_HC_accum(d_val,h_val,age_j,xi_1,xi_2,del_h,h_l);
        [ind_l,weight_l] = find_loc(h_grid,hprime_val);
        ind_r    = ind_l+1;
        weight_r = 1-weight_l;
        for zp_c=1:n_z
            mu(aprime,ind_l,zp_c,j) = mu(aprime,ind_l,zp_c,j)+weight_l*pi_z(z_c,zp_c)*mu(a_c,h_c,z_c,j-1);
            mu(aprime,ind_r,zp_c,j) = mu(aprime,ind_r,zp_c,j)+weight_r*pi_z(z_c,zp_c)*mu(a_c,h_c,z_c,j-1);
        end
        %mu(aprime,ind_l,:,j) = mu(aprime,ind_l,:,j)+shiftdim(weight_l*pi_z(z_c,:)*mu(a_c,h_c,z_c,j-1),-1);
        %mu(aprime,ind_r,:,j) = mu(aprime,ind_r,:,j)+shiftdim(weight_r*pi_z(z_c,:)*mu(a_c,h_c,z_c,j-1),-1);
    end
    end
    end

    if verbose == 1
        fprintf('Age = %d out of %d \n',j,N_j)
        check = sum(mu(:,:,:,j),"all");
        disp(check)
    end

end %end j

%% Multiply by age distribution
mu_j = Params.mewj;

for j=1:N_j
    mu(:,:,:,j) = mu(:,:,:,j)*mu_j(j);
end

StationaryDist = mu;

end %end function