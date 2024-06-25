classdef fun
    %This class contains all functions related to the model (functional
    %forms)

    methods (Static)


        % FnsToEvaluate.l_coh=@(l_f,aprime,a,h_f,eta_m,eta_f,theta) l_f;
        % FnsToEvaluate.a_coh=@(l_f,aprime,a,h_f,eta_m,eta_f,theta) a;
        % FnsToEvaluate.h_coh=@(l_f,aprime,a,h_f,eta_m,eta_f,theta) exp(h_f);
        % FnsToEvaluate.ym_coh=@(l_f,aprime,a,h_f,eta_m,eta_f,theta,w_m,eff_j) w_m*eff_j*theta*eta_m;
        % FnsToEvaluate.yf_coh=@(l_f,aprime,a,h_f,eta_m,eta_f,theta,w_f) w_f*exp(h_f)*theta*eta_f*l_f;
        % FnsToEvaluate.pen_coh=@(l_f,aprime,a,h_f,eta_m,eta_f,theta,pen_j) pen_j;
        % FnsToEvaluate.y_coh = @(l_f,aprime,a,h_f,eta_m,eta_f,theta,w_m,eff_j,w_f) w_m*eff_j*theta*eta_m+w_f*exp(h_f)*theta*eta_f*l_f;
        % FnsToEvaluate.c_coh=@(l_f,aprime,a,h_f,eta_m,eta_f,theta,w_m,eff_j,w_f,pchild_j,pen_j,r,agej,Jr) f_ReturnFn_cons(l_f,aprime,a,h_f,eta_m,eta_f,theta,w_m,eff_j,w_f,pchild_j,pen_j,r,agej,Jr);
        % FnsToEvaluate.c_eq_coh = @(l_f,aprime,a,h_f,eta_m,eta_f,theta,w_m,eff_j,w_f,pchild_j,pen_j,r,nchild_j,agej,Jr) f_ReturnFn_cons_equiv(l_f,aprime,a,h_f,eta_m,eta_f,theta,w_m,eff_j,w_f,pchild_j,pen_j,r,nchild_j,agej,Jr);
        %


        function F = l_coh(l_f,aprime,a,h_f,eta_m,eta_f,theta)
            F = l_f;
        end

        function F = a_coh(l_f,aprime,a,h_f,eta_m,eta_f,theta)
            F = a;
        end

        function F = h_coh(l_f,aprime,a,h_f,eta_m,eta_f,theta)
            F = exp(h_f);
        end

        function F = ym_coh(l_f,aprime,a,h_f,eta_m,eta_f,theta,w_m,eff_j)
            F = w_m*eff_j*theta*eta_m;
        end

        function F = yf_coh(l_f,aprime,a,h_f,eta_m,eta_f,theta,w_f)
            F = w_f*exp(h_f)*theta*eta_f*l_f;
        end

        function F = y_coh(l_f,aprime,a,h_f,eta_m,eta_f,theta,w_m,eff_j,w_f)
            F = w_m*eff_j*theta*eta_m+w_f*exp(h_f)*theta*eta_f*l_f;
        end

        function F = pen_coh(l_f,aprime,a,h_f,eta_m,eta_f,theta,pen_j)
            F = pen_j;
        end

        function F = c_coh(l_f,aprime,a,h_f,eta_m,eta_f,theta,w_m,eff_j,w_f,pchild_j,pen_j,r,agej,Jr)
            F = f_ReturnFn_cons(l_f,aprime,a,h_f,eta_m,eta_f,theta,w_m,eff_j,w_f,pchild_j,pen_j,r,agej,Jr);
        end

        function F = c_eq_coh(l_f,aprime,a,h_f,eta_m,eta_f,theta,w_m,eff_j,w_f,pchild_j,pen_j,r,nchild_j,agej,Jr)
            F = f_ReturnFn_cons_equiv(l_f,aprime,a,h_f,eta_m,eta_f,theta,w_m,eff_j,w_f,pchild_j,pen_j,r,nchild_j,agej,Jr);
        end


    end %END METHODS
end %END CLASS <fun>