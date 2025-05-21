freq_R=[35.66, 37.54, 39.37, 42.81, 46.99, 33.03, 31.25, 28.12, 24.13, 20.34, 100.7, 10.02]*1e3;
R=[10,100,220,470, 1000];
err_R=[13.8,99.6,220.3,465,1004];

ampl_diff = zeros(size(file_names));
phase_diff = zeros(size(file_names));
err_amp = zeros(size(file_names));
err_phase = zeros(size(file_names));
for i = 1:size(file_names, 2)
    data1 = readmatrix(strcat('470ohm-b/', file_names(i), '/meas1.csv'));
    data2 = readmatrix(strcat('470ohm-b/', file_names(i), '/meas2.csv'));

    % Model function to fit: y = a * exp(b * x)
    if isempty(find(indexes1k == i))
        freq = freq100k(find(indexes100k == i));
    else
        freq = freq1k(find(indexes1k == i));
    end
    model = @(params, x) params(1) * cos(2*pi*freq * x + params(2));

    lb = [0.3, -pi];
    ub = [0.5, pi];

    [paramsFitted1, ~, residual, ~, ~, ~, jacobian] = lsqcurvefit(model, [0.4, 0], data1(:, 1), data1(:, 2), lb, ub);
    var_res = sum(residual.^2) / (length(residual) - length(paramsFitted1));
    cov_matrix = var_res * inv(jacobian' * jacobian);
    std_errors_1 = full(sqrt(diag(cov_matrix)));
    [paramsFitted2, ~, residual, ~, ~, ~, jacobian] = lsqcurvefit(model, [0.02, 0], data2(:, 1), data2(:, 2), [0, -pi], ub);
    var_res = sum(residual.^2) / (length(residual) - length(paramsFitted2));
    cov_matrix = var_res * inv(jacobian' * jacobian);
    std_errors_2 = full(sqrt(diag(cov_matrix)));

    ampl_diff(i) = abs(paramsFitted2(1) / paramsFitted1(1));
    phase_diff(i) = mod(paramsFitted1(2) - paramsFitted2(2), pi);
    err_amp(i) = scope_err(data1(1, 3)) + scope_err(data2(1, 3));
    err_phase(i) = std_errors_1(2) + std_errors_2(2);
end

clear alpha beta cov_matrix delta dx Fs gamma i jacobian k L lb model N P paramsFitted1 paramsFitted2 refinedFreq reifnedIndex residual std_errors_1 std_errors_2 ub var_res w xdata Y y_detrended y_windowed ydata B0 data1 data2
