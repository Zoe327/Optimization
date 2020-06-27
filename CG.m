% Conjugate Gradient
function [ V_min, X_opt, iter, v_set, xj_set ] = CG(f,v,x0, flag, threshold)
% initialization
j = 0;
xj = x0;
v_set = [];   % v_set and xj_set are writted to plot convergence curve
xj_set = xj;
con = true;
% calculate gradient  and define search direction
% with jacobian function or finite difference algorithm
if (flag == 0)
    grad_V = jacobian(f,v);
    % define the search direction
    sj = -subs(grad_V, v, xj);
else
    sj = -e_finiteDiff(f,v,xj);
end

% since iteration 1, introduce beta to update search direction
while (con)
    
    % stop condition
    if (norm(sj) < threshold)
        con = false;
    else
        % find the optimal step size using Armijo algorithm
        w = b_Armijo(f,v,sj',xj, 1.5, 0.8);
        % update the point
        v_set = [v_set subs(f,v,xj)];
        xi_old = xj;
        xj = xi_old + w*sj';
        double(xj)
        xj_set = [xj_set xj];
        
        % calculate the gradient according to the method specified
        if (flag == 0)
            grad_new = subs(grad_V,v,xj);
            grad_old = subs(grad_V,v,xi_old);
            
        else
            grad_new = e_finiteDiff(f,v,xj);
            grad_old = e_finiteDiff(f,v,xi_old);
        end
        
        beta = ((grad_new)-(grad_old))*(grad_new)'/...
            ((grad_old)*(grad_old)') ;     
        
        sj_old = sj;
        sj = -grad_new + beta*sj_old;
        j = j + 1;
    end
    
end

V_min = subs(f,v,xj);
X_opt = xj;
iter = j;
end