% Steepest Descent 
function [ v_min, x_opt, iter, V_set, x_set ]...
    = SD(f,v,x0, flag, threshold)
% initialization 
k = 0;
xj = x0;
V_set = [];
x_set = xj;
con = true;
% flag = 0 use jacobian to calculate gradient
% flag = 1 use finite difference algorithm to claculate gradient
if (flag == 0)
    grad_V = jacobian(f,v);
end 

while (con)
    if (flag == 1) 
        sj = -finiteDiff(f,v,xj)
        V_set = [V_set subs(f,v,xj)];
    else 
        sj = -subs(grad_V, v, xj)   % search direction 
        V_set = [V_set subs(f,v,xj)];
    end
    
    % stop condition 
    if (norm(sj) < threshold)
        con = false;
    else  % update point
        w = Armijo(f,v,sj',xj, 1.5, 0.8);  % Armijo Agorithm - stepsize
        xj = xj + w*sj';  
        x_set = [x_set xj];
        k = k + 1;
    end
end
v_min = subs(f,v,xj);
x_opt = xj;
iter = k;
end
