function [z,pi] = discretize_AR(rho,mu,sigma_eps,n)

% calculate variance of the overall process
sigma_eta = sigma_eps/(1-rho^2);

% determine the transition matrix
pi = rouwenhorst_matrix(rho,n);

%determine the nodes
psi = sqrt(n-1)*sqrt(sigma_eta);
z = zeros(n,1);
for in = 1:n
    z(in) = -psi + 2*psi*(in-1)/(n-1);
end
z = z + mu;

end %end function

function [pi_new] = rouwenhorst_matrix(rho,n)


p = (1 + rho)/2;

if (n==2)
    pi_new = zeros(n,n);
    pi_new(1, :) = [p, 1-p];
    pi_new(2, :) = [1-p, p];
else
    pi_old = rouwenhorst_matrix(rho,n);
    pi_new = zeros(n,n);

    pi_new(1:n-1, 1:n-1) = pi_new(1:n-1, 1:n-1) + p*pi_old;
    pi_new(1:n-1, 2:n  ) = pi_new(1:n-1, 2:n  ) + (1-p)*pi_old;
    pi_new(2:n  , 1:n-1) = pi_new(2:n  , 1:n-1) + (1-p)*pi_old;
    pi_new(2:n  , 2:n  ) = pi_new(2:n  , 2:n  ) + p*pi_old;

    pi_new(2:n-1, :) = pi_new(2:n-1, :)/2;

end

end %end function rouwenhorst_matrix

%     recursive subroutine rouwenhorst_matrix(rho, pi_new)
% 
%         implicit none
%         real(8), intent(in) :: rho
%         real(8), intent(out) :: pi_new(:, :)
%         integer :: n
%         real(8) :: p, pi_old(size(pi_new,1)-1, size(pi_new,1)-1)
% 
%         n = size(pi_new, 1)
%         p = (1d0 + rho)/2d0
% 
%         if(n == 2)then
%             pi_new(1, :) = [p, 1d0-p]
%             pi_new(2, :) = [1d0-p, p]
%         else
%             call rouwenhorst_matrix(rho, pi_old)
%             pi_new = 0d0
% 
%             pi_new(1:n-1, 1:n-1) = pi_new(1:n-1, 1:n-1) + p*pi_old
%             pi_new(1:n-1, 2:n  ) = pi_new(1:n-1, 2:n  ) + (1d0-p)*pi_old
%             pi_new(2:n  , 1:n-1) = pi_new(2:n  , 1:n-1) + (1d0-p)*pi_old
%             pi_new(2:n  , 2:n  ) = pi_new(2:n  , 2:n  ) + p*pi_old
% 
%             pi_new(2:n-1, :) = pi_new(2:n-1, :)/2d0
%         endif
%     end subroutine
% 


% !##############################################################################
%     ! SUBROUTINE discretize_AR
%     !
%     ! Discretizes an AR(1) process of the form z_j = \rho*z_{j-1} + eps using
%     !     the Rouwenhorst method.
%     !
%     ! REFERENCE: Kopecky, K.A., Suen, R.M.H., Finite state Markov-chain 
%     !            approximations to highly persistent processes, Review of Economic
%     !            Dynamics, Vol. 13, No. 3, 2010, 701-714.
%     !##############################################################################
% subroutine discretize_AR(rho, mu, sigma_eps, z, pi, w)
% 
%     implicit none
%     ! ------------------------- DECLARE INPUTS ---------------------------------!
%     ! autoregression parameter
%     real(8), intent(in) :: rho
%     ! unconditional mean of the process
%     real(8), intent(in) :: mu
%     ! variance (NOT standard deviation!) of the shock
%     real(8), intent(in) :: sigma_eps
%     ! number of points
%     ! This is size(z)
%     ! ------------------------- DECLARE OUTPUTS ---------------------------------!
%     ! discrete shock values
%     real(8), intent(out) :: z(:)
%     ! transition matrix
%     real(8), intent(out) :: pi(:, :)
%     ! the stationary distribution
%     real(8), intent(out), optional :: w(:)
%     ! ------------------------- DECLARE LOCALS ---------------------------------!
%     integer :: n, in
%     real(8) :: psi, sigma_eta
% 
%     ! assert size equality and get approximation points
%     n = size(z)
%     if (n /= size(pi,1)) then
%         call myerror("discretize_AR: z and pi do not match")
%     endif
%     if (size(pi,1) /= size(pi,2)) then
%         call myerror("discretize_AR: pi is not square matrix")
%     endif
% 
%     ! calculate variance of the overall process
%     sigma_eta = sigma_eps/(1d0-rho**2)
% 
%     ! determine the transition matrix
%     call rouwenhorst_matrix(rho, pi)
% 
%     ! determine the nodes
%     psi = sqrt(real(n-1,8))*sqrt(sigma_eta)
%     do in = 1,n
%         z(in) = -psi + 2d0*psi*real(in-1,8)/real(n-1,8)
%     enddo
%     z = z + mu
% 
%     if (present(w)) then
%         w = 1d0/real(n,8)
%         do in = 1, 10000
%             w = matmul(transpose(pi), w)
%         enddo
%     endif
% 
%     !##########################################################################
%     ! Subroutines and functions
%     !##########################################################################
% 
%     contains
% 
%     !##########################################################################
%     ! subroutine rouwenhorst_matrix
%     !
%     ! Calculates value of function that should be integrated for pis.
%     !##########################################################################
%     recursive subroutine rouwenhorst_matrix(rho, pi_new)
% 
%         implicit none
%         real(8), intent(in) :: rho
%         real(8), intent(out) :: pi_new(:, :)
%         integer :: n
%         real(8) :: p, pi_old(size(pi_new,1)-1, size(pi_new,1)-1)
% 
%         n = size(pi_new, 1)
%         p = (1d0 + rho)/2d0
% 
%         if(n == 2)then
%             pi_new(1, :) = [p, 1d0-p]
%             pi_new(2, :) = [1d0-p, p]
%         else
%             call rouwenhorst_matrix(rho, pi_old)
%             pi_new = 0d0
% 
%             pi_new(1:n-1, 1:n-1) = pi_new(1:n-1, 1:n-1) + p*pi_old
%             pi_new(1:n-1, 2:n  ) = pi_new(1:n-1, 2:n  ) + (1d0-p)*pi_old
%             pi_new(2:n  , 1:n-1) = pi_new(2:n  , 1:n-1) + (1d0-p)*pi_old
%             pi_new(2:n  , 2:n  ) = pi_new(2:n  , 2:n  ) + p*pi_old
% 
%             pi_new(2:n-1, :) = pi_new(2:n-1, :)/2d0
%         endif
%     end subroutine
% 
% end subroutine discretize_AR
% !===============================================================================!
