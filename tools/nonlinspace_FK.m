function a_grid = nonlinspace_FK(a_min,a_max,n,growth)
    % Fehr and Kindermann
    % Build a non-equally spaced grid 
    % a_min:  lower bound
    % a_max:  upper bound
    % n:      number of points
    % growth: parameter to govern the spacing of grid points
    
   if n>=2
        % Kindermann's method
        a_grid = nan(n,1);
        % calculate factor
        h      = (a_max-a_min)/((1+growth)^n-1);
        for i=2:n+1
            a_grid(i-1) = h*((1+growth)^(i-1)-1) + a_min;
        end
    else
        a_grid = a_min; % as in Matlab linspace
    end
    
end