function [V,Policy] = ValueFnIter_Case1_FHorz_ale(n_d,n_a,n_z,N_j, d_grid, a_grid, z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames,XXX, vfoptions,par)
%% THIS RUNS ON CPU
% STATE VARIABLES
% V(a,h,z,j)
% a: asset holdings
% h: human capital (female)
% z: labor productivity shocks (eta_m, eta_f are shocks, theta is permanent)
% j: age (from 1 to N_j)
% CHOICE VARIABLES
% d:  Female labor supply (only extensive margin: either 0 or 1)
% a': Next-period assets
% h' is implied by (d,a) and consumption is implied by the budget
% constraint
% DYNAMIC PROGRAMMING PROBLEM
% V(a,h,z,j) = max_{d,a'} F(d,a',a,h,z,j)+beta*s_j*E[V(a',h',z',j+1)|z]
% subject to
% h'=G(d,h), law of motion for human capital

verbose = vfoptions.verbose;
do_refine = par.do_refine;
n_ref     = par.n_ref;

%% Define grids and grid sizes
[n_a,n_h,a_grid,h_grid,z_grid,n_z] = prepare_grids(n_a,a_grid,z_grid,n_z,pi_z);

aprime_val = a_grid;  %(a',1), values
a_val      = a_grid'; %(1,a), values
a_grid_row = a_grid'; %(1,a), values
a_index    = (1:n_a)'; %(a,1), index

%% Set parameters that do not depend on age
beta = Params.beta;
r    = Params.r;
w_m  = Params.w_m;
w_f  = Params.w_f;
crra = Params.crra;
nu   = Params.nu;
Jr   = Params.Jr;
xi_1  = Params.xi_1;
xi_2  = Params.xi_2;
del_h = Params.del_h;
h_l   = Params.h_l;

% Initialize output arrays
V  = zeros(n_a,n_h,n_z,N_j);
Policy = zeros(2,n_a,n_h,n_z,N_j);

%% Solve problem in the last period

% Set age-dependent parameters
eff_j     = Params.eff_j(N_j);
pchild_j  = Params.pchild_j(N_j);
pen_j     = Params.pen_j(N_j);
nchild_j  = Params.nchild_j(N_j);
age_j     = Params.agej(N_j);

% V(a,h,z,N_j) = max_{d,a'}
V_d = zeros(n_a,n_h,n_z,n_d);
Pol_aprime_d = zeros(n_a,n_h,n_z,n_d);

for d_c = 1:n_d
    d_val = d_grid(d_c);
    for z_c = 1:n_z
        eta_m_val = z_grid(z_c,1);
        eta_f_val = z_grid(z_c,2);
        theta_val = z_grid(z_c,3);
        for h_c = 1:n_h
            h_val = h_grid(h_c);
            % RetMat is (a',a)
            RetMat = ReturnFn_ale(d_val,aprime_val,a_val,h_val,eta_m_val,eta_f_val,theta_val,...
                w_m,eff_j,w_f,pchild_j,pen_j,r,nchild_j,crra,nu,age_j,Jr);
            [max_val,max_ind] = max(RetMat,[],1);
            V_d(:,h_c,z_c,d_c) = max_val;
            Pol_aprime_d(:,h_c,z_c,d_c) = max_ind;
        end %end h
    end %end z
end %end d

[V(:,:,:,N_j),d_max] = max(V_d,[],4);

Policy(1,:,:,:,N_j) = d_max; % Optimal d
for z_c=1:n_z
    for h_c=1:n_h
        for a_c = 1:n_a
            d_star = d_max(a_c,h_c,z_c);
            Policy(2,a_c,h_c,z_c,N_j) = Pol_aprime_d(a_c,h_c,z_c,d_star); % Optimal a'
        end
    end
end

%% Backward iteration over age

for j = N_j-1:-1:1

    if verbose==1; fprintf('Age %d out of %d \n',j,N_j); end

    V_next = V(:,:,:,j+1); %V(a',h',z')

    % Set age-dependent parameters
    eff_j     = Params.eff_j(j);
    pchild_j  = Params.pchild_j(j);
    pen_j     = Params.pen_j(j);
    nchild_j  = Params.nchild_j(j);
    age_j     = Params.agej(j);
    s_j       = Params.s_j(j);

    for z_c = 1:n_z

        eta_m_val = z_grid(z_c,1);
        eta_f_val = z_grid(z_c,2);
        theta_val = z_grid(z_c,3);

        % Compute EV(a',h'), given z
        % EV = zeros(n_a,n_h);
        z_prob = pi_z(z_c,:)';
        EV = V_next.*shiftdim(z_prob,-2); %V(a',h',z')*Prob(1,1,z')
        EV = sum(EV,3); %EV(a',h')

        for d_c=1:n_d
            d_val = d_grid(d_c);
            for h_c = 1:n_h
                h_val = h_grid(h_c);
                hprime_val = f_HC_accum(d_val,h_val,age_j,xi_1,xi_2,del_h,h_l);
                % Ret_mat is (a',a) since aprime_val is (a',1) and a_val is
                % (1,a)
                Ret_mat = ReturnFn_ale(d_val,aprime_val,a_val,h_val,eta_m_val,eta_f_val,theta_val,...
                    w_m,eff_j,w_f,pchild_j,pen_j,r,nchild_j,crra,nu,age_j,Jr);
                %[ind_l,weight_l] = interp_toolkit_scal(hprime_val,h_grid);
                [ind_l,weight_l] = find_loc(h_grid,hprime_val);
                EV_interp = EV(:,ind_l)*weight_l+EV(:,ind_l+1)*(1-weight_l);
                RHS_mat = Ret_mat+beta*s_j*EV_interp;
                [max_val,max_ind] = max(RHS_mat,[],1);
                % max_val and max_ind are (1,a)
                V_d(:,h_c,z_c,d_c) = max_val;  %best V given d
                Pol_aprime_d(:,h_c,z_c,d_c) = max_ind; %best a' given d, index

                if do_refine==1
                    Interpolant = griddedInterpolant(a_grid,EV_interp);
                    ind_l = max(max_ind-1,1);   %(1,a)
                    ind_r = min(max_ind+1,n_a); %(1,a)
                    % a_grid_fine is (a'_ref,n_a)
                    a_grid_fine = nonlinspace(a_grid_row(ind_l),a_grid_row(ind_r),n_ref,1.0);
                    % Ret_mat_refine is (a'_ref,n_a)
                    % since a_grid_fine is (a'_ref,a) and a_val is (1,a)
                    % d_val and h_val are scalars
                    Ret_mat_refine = ReturnFn_ale(d_val,a_grid_fine,a_val,h_val,eta_m_val,eta_f_val,theta_val,...
                        w_m,eff_j,w_f,pchild_j,pen_j,r,nchild_j,crra,nu,age_j,Jr);
                    % RHS_mat_refine is (a'_ref,a)
                    RHS_mat_refine = Ret_mat_refine+beta*s_j*Interpolant(a_grid_fine);
                    % max_val and max_ind are (1,a)
                    [max_val,max_ind] = max(RHS_mat_refine,[],1);
                    V_d(:,h_c,z_c,d_c) = max_val;  %best V given d
                    % TODO: max_ind now refers to a_grid_refine, not
                    % a_grid! Need interpolation over a' when doing the
                    % distribution
                    Pol_aprime_d(:,h_c,z_c,d_c) = max_ind; %best a' given d, index
                end %end if    
            end %end h
        end %end d
    end %end z

    [V(:,:,:,j),d_max] = max(V_d,[],4);

    Policy(1,:,:,:,j) = d_max; % Optimal d
    for z_c=1:n_z
        for h_c=1:n_h
            for a_c = 1:n_a
                d_star = d_max(a_c,h_c,z_c);
                Policy(2,a_c,h_c,z_c,j) = Pol_aprime_d(a_c,h_c,z_c,d_star); % Optimal a'
            end
        end
    end

end %end j

% % Transform 4-dim Policy(a,h,z,j) in 6-dim Policy_out(a,h,eta_m,eta_f,theta,j)
% % Policy_out = zeros(2,n_a,n_h,length(eta_m_grid),length(eta_f_grid),length(theta_grid),N_j);
% % for j=1:N_j
% %     for z_c=1:n_z
% %         for h_c=1:n_h
% %             for a_c = 1:n_a
% %                 [eta_m_c,eta_f_c,theta_c]=ind2sub([length(eta_m_grid),length(eta_f_grid),length(theta_grid)],z_c);
% %                 Policy_out(1,a_c,h_c,eta_m_c,eta_f_c,theta_c,j)=Policy(1,a_c,h_c,z_c,j);
% %                 Policy_out(2,a_c,h_c,eta_m_c,eta_f_c,theta_c,j)=Policy(2,a_c,h_c,z_c,j);
% %             end
% %         end
% %     end
% % end
% Policy_out = reshape(Policy,[2,n_a,n_h,length(eta_m_grid),length(eta_f_grid),length(theta_grid),N_j]);

end %end function