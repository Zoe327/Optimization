% Secant Algorithm
function [ V_opt, x_opt, Iter, v_set, x_set] = Secant( f,v,x0, flag, threshold )

wthreshold = 10^-4;
H = eye(length(v));  %identical matrix
j = 0;
xj = x0;
run = true;
v_set = [];
x_set = xj;
%n=0;

if (flag == 0)
    grad_V = jacobian(f,v);  % define gradient 
end

while (run)
    % find search direction 
    if (flag == 0)
        if (j == 0)
            grad_V_value = subs(grad_V, v, xj);
        else
            grad_V_value = grad_xj;
        end
        
        sj = -H*(grad_V_value)';
    else
        sj = -H*(e_finiteDiff(f,v,xj))';
    end
    
    v_set = [v_set subs(f,v,xj)];
    
    w = b_Armijo(f,v,sj,xj,1.5,0.8);  % step size 
    % update point and H with DFP
    old_xj = xj;
    xj = old_xj + w*sj;
    x_set = [x_set xj];
    if (flag == 0)
        grad_xj = subs(grad_V, v, xj);
        grad_old_xj = grad_V_value;
    else
        grad_xj = e_finiteDiff(f,v,xj);
        grad_old_xj = e_finiteDiff(f,v,old_xj);
    end
    % check if the gradient is smaller than the threshold
    % if yes, then exit the loop, else, update the H matrix
    if (double(norm(grad_xj)) < threshold)
        run = false;
    else
        delta_xj = xj - old_xj;
        delta_gi = grad_xj' - grad_old_xj';
        
        H = H + ((delta_xj*delta_xj')./(delta_xj'*delta_gi))...
            - ((H*delta_gi)*(H*delta_gi)')./(delta_gi'*H*delta_gi);
        j = j+1;
    end
    a = double(subs(f,v,xj));
    b = double(subs(f,v,old_xj));
    if (a>b)
        H = eye(length(v));
        fprintf('restart');
    end
  %  n=n+1
end
V_opt = subs(f,v,xj);
x_opt = xj;
Iter = j;
end
