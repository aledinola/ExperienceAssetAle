function [n_a,n_h,a_grid,h_grid,z_grid,n_z] = prepare_grids(n_a,a_grid,z_grid,n_z,pi_z)

n_a_all    = n_a;
a_grid_all = a_grid;
z_grid_all = z_grid;
n_a = n_a_all(1);
n_h = n_a_all(2);
a_grid = a_grid_all(1:n_a_all(1));
h_grid = a_grid_all(n_a_all(1)+1:n_a_all(1)+n_a_all(2));
eta_m_grid = z_grid_all(1:n_z(1));
eta_f_grid = z_grid_all(n_z(1)+1:n_z(1)+n_z(2));
theta_grid = z_grid_all(n_z(1)+n_z(2)+1:end);

[x1,x2,x3] = ndgrid(eta_m_grid,eta_f_grid,theta_grid);
x1 = x1(:);
x2 = x2(:);
x3 = x3(:);
z_grid = [x1,x2,x3];
n_z = size(z_grid,1);
if n_z~=size(pi_z,1) || n_z~=size(pi_z,2)
    error('n_z and pi_z are not consistent')
end

end %end function