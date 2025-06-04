function [params_fit] = do_fit(fun, x, y, startpoint)
    lb = zeros(size(starpoint));
    ub = ones(size(startpoint)) * Inf;

    loss_fun = @(p) loss_function(y, fun(p, x));
    options = optimoptions('fmincon', 'Display', 'iter');
    [params_fit, ~] = fmincon(loss_fun, startpoint, [], [], [], [], lb, ub, [], options);
end

function [loss] = loss_function(real_y, fitted_y)
    loss = sum((real_y - fitted_y).^2);
end