function [ V_min, x_opt, Iter, Cxeq ] = ...
    LagNewton( f,v, x0, Ceq, threshold, lambda )

% initialize the parameteres
xj = x0;
run = true;
j = 1;
ithreshold = 10^-4;
xj_set = x0;

while run
    % construct the Lagrangian
    L = f - lambda.'*Ceq;
    
    % calculate the elements required for the matrix operation
    gradLx = jacobian(L,v);
    gradLxx = subs(jacobian(gradLx,v), v, xj);    
    gradCx = subs(jacobian(Ceq, v), v, xj);
    gradfx = subs(jacobian(f,v), v, xj);  
    CeqValue = subs(Ceq, v, xj)
    
    % find the updates
    deltaX = -1.*(gradCx)\CeqValue
    deltaLambda = gradCx'\(-1*gradLxx*deltaX + gradfx');
    
    if (norm(CeqValue) < threshold)
        run = false;
    else
        % update the parameters accordingly
        lambda = lambda + deltaLambda;
        xj = xj + deltaX;
        xj_set = [xj_set xj];
        j = j+1;
    end
    
end

V_min = subs(f,v,xj);
x_opt = xj;
Iter = j;
end
