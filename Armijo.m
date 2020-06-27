
function [stepsize] = Armijo( f, v, sj, xj, r, u)

% initialize the power of gamma and mu
p = 0;
q = 0;

if (r <= 1)
    r = 1.5;
end
if (u >= 1 || u <= 0)
    u = 0.8;
end

% compute the gradient of the cost function f wrt variables v
grad_f = jacobian(f,v);
w = r^p;
con = true;

% Find a w on the right of wj
while (con)
    % calculate the functions
    f_bar = subs(f,v,xj)+(1/2)*w*(subs(grad_f,v,xj)*sj);
    step = xj + w*sj;
    f_value = subs(f,v,step);
%     k = f_bar - f_value;
    % if the condition is met, stop, otherwise, increment the order p
    if (f_value >= f_bar)
        con = false;
    else
        p = p + 1;
        w = r^p;
    end
end

% re-calculate w
w = (r^p)*(u^q);

% Find a w on the left of wj
while (not(con))
    % calculate the functions
    f_bar = subs(f,v,xj)+(1/2)*w*(subs(grad_f,v,xj)*sj);
    step = xj + w*sj;
    f_value = subs(f,v,step);
    %     k = f_value - f_bar
    % if the condition is met, stop, otherwise, increment the order q
    if (f_value <= f_bar)
        con = true;
    else
        q = q+1;
        w = (r^p)*(u^q);
    end
end
stepsize = w;
end
