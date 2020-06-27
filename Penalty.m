% penalty 
function [ V_min, x_opt, Iter, cx_set, gx_set, xj_set ]...
    = Penalty( V,x, Ceq, Cineq, sigma, x0 , threshold)
% initilization
j = 1;
xj = x0;
CON1 = true;   % whether find to optimal value
CON2 = true;   % whether met the inequality constraints
n = length(Cineq);
xj_set = xj;
cx_set = [];
gx_set = [];

while CON1
    
    sigma = sigma^j;  % penalty
    
    % calculate the value of ineq at xj
    gx_i = subs(Cineq, x, xj);
    gx_set = [gx_set gx_i];
    ModCineq = Cineq;
    % if the g(x) <= 0, then the ineq is met, remove the corresponding constraint
    for i = 1:n
        if (gx_i(i,1) <= 0)
            ModCineq(i,1) = 0;
        end
    end
    
    % then the constraint problem can be converted to an unconstraint one
    % solve the unconstrainted problem using secant Algorithm
    Phi = V + sigma*(Ceq.^2 + ModCineq.^2);
    [Vi_min, Xi_opt, Iter, vi_set, xji_set] =...
        Secant_debug(Phi, x, xj, 1, 10^-3);
    xj = Xi_opt
    xj = double(xj); 
    xj_set = [xj_set Xi_opt];
   
    % case where there is no equality constraint is automatically
    % handled by having a sum smaller than threshold
    c_j = subs(Ceq,x,xj);
    cx_set = [cx_set c_j];
    if (sum(abs(c_j)) < threshold)
        g_x = subs(Cineq,x,xj);
        % if the equality condition is met (by having a value smaller
        % than threshold, check the inequality condition
        for i=1:n
            if (g_x(i,1) > 0)
                CON2 = false;
            end
            CON2 = true;
        end
        % if the inequality condition is met, then terminate
        if (CON2)
            CON1 = false;
        end
        
    else
        j = j+1
    end
end
x_opt = xj;
Iter = j;
V_min = subs(V,x,xj);
end