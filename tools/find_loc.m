function [jl,omega] = find_loc(x_grid,xi)

%-------------------------------------------------------------------------%
% DESCRIPTION
% Find jl s.t. x_grid(jl)<=xi<x_grid(jl+1)
% for jl=1,..,N-1
% omega is the weight on x_grid(jl) so that
% omega*x_grid(jl)+(1-omega)*x_grid(jl+1)=xi
% INPUTS
% x_grid must be a strictly increasing column vector (nx,1)
% xi must be a scalar
% OUTPUTS
% jl: Left point (scalar)
% omega: weight on the left point (scalar)
% NOTES
% See find_loc_vec.m for a vectorized version.
%-------------------------------------------------------------------------%

nx = size(x_grid,1);

jl = max(min(locate(x_grid,xi),nx-1),1);
%Weight on x_grid(j)
omega = (x_grid(jl+1)-xi)/(x_grid(jl+1)-x_grid(jl));
omega = max(min(omega,1),0);

end %end function "find_loc"