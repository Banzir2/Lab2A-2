close all; clear;
freq_H=[33.8,36.39, 38.58, 41.25, 45.17, 31.13, 29.38, 27.29, 23.19, 105.7, 10.15]*1e3;
file_names = ["34k", "36k", "38k", "41k", "45k", "31k", "27k", "29k", "23k","106k","10k"];

w = 2 * pi * freq_H(:);    % Angular frequencies

% Initialize
ampl_diff = zeros(size(file_names));
phase_diff = zeros(size(file_names));
err_amp = zeros(size(file_names));
err_phase = zeros(size(file_names));

for i = 1:length(file_names)
    data1 = readmatrix(strcat('470ohm-c/', file_names(i), '/meas1.csv'));
    data2 = readmatrix(strcat('470ohm-c/', file_names(i), '/meas2.csv'));

    % Model function to fit: y = a * exp(b * x)
    freq = freq_H(i);
    model = @(params, x) params(1) * cos(2*pi*freq * x + params(2));

    lb = [0, -pi];
    ub = [25, pi];

    [paramsFitted1, ~, residual, ~, ~, ~, jacobian] = lsqcurvefit(model, [0.4, 0], data1(:, 1), data1(:, 2), lb, ub);
    var_res = sum(residual.^2) / (length(residual) - length(paramsFitted1));
    cov_matrix = var_res * inv(jacobian' * jacobian);
    std_errors_1 = full(sqrt(diag(cov_matrix)));
    [paramsFitted2, ~, residual, ~, ~, ~, jacobian] = lsqcurvefit(model, [0.02, 0], data2(:, 1), data2(:, 2), [0, -pi], ub);
    var_res = sum(residual.^2) / (length(residual) - length(paramsFitted2));
    cov_matrix = var_res * inv(jacobian' * jacobian);
    std_errors_2 = full(sqrt(diag(cov_matrix)));

    ampl_diff(i) = abs(paramsFitted2(1) / paramsFitted1(1));
    raw_diff = paramsFitted1(2) - paramsFitted2(2);
    phase_diff(i) = raw_diff;    err_amp(i) = scope_err(data1(1, 3)) + scope_err(data2(1, 3));
    err_phase(i) = std_errors_1(2) + std_errors_2(2);
end
phase_diff = unwrap(phase_diff);
%Amp diff
function r = weighted_residuals(model, x, xdata, ydata, weights)
    r = weights .* (model(x, xdata) - ydata);
end
xhisquare = zeros(1, 2);
C = 1.0571e-09;
H = 21.06e-3;

% --- Model and Fit ---
modelFun = @(params, w) (params(1) .* w) ./ sqrt(1 + (params(1).*w - params(2)./w).^2);
startPoint = [H/470, 1/(470*C)];         % [L/R, 1/RC]
xData = 2*pi*freq_H(:);
yData = ampl_diff(:);
lb = [0, 0];
ub = [Inf, Inf];
options = optimoptions('lsqcurvefit', 'TolFun', 1e-8, 'Display', 'off');

% Fit model
[params_fit, ~, residuals, ~, ~, ~, jacobian] = lsqcurvefit(modelFun, startPoint, xData, yData, lb, ub, options);

% --- Error Estimation (Confidence Intervals) ---
var_res = sum(residual.^2) / (length(residual) - length(params_fit));
cov_matrix = var_res * inv(jacobian' * jacobian);
param_errors = full(sqrt(diag(cov_matrix)));

% --- Evaluate model over a fine frequency range ---
freq_plot = linspace(min(freq_H), max(freq_H), 500)';
x_plot = 2 * pi * freq_plot;
y_fit = modelFun(params_fit, x_plot);

% --- Error propagation for the model ---
% Approximate uncertainty band (using finite differences)
delta_p1 = param_errors(1);
delta_p2 = param_errors(2);
y_high = modelFun(params_fit + [delta_p1, delta_p2], x_plot);
y_low  = modelFun(params_fit - [delta_p1, delta_p2], x_plot);
y_err = abs(y_high - y_low) / 2;   % Symmetric error estimate

% --- Reduced Chi-Squared ---
N = length(yData);        % number of data points
k = length(params_fit);   % number of parameters
reduced_chi2 = sum((residuals ./ err_amp(:)).^2) / (N - k);

% --- Plotting ---
figure; hold on;

% Plot data with error bars
errorbar(freq_H, yData, err_amp, '.', 'MarkerSize', 25, 'DisplayName', 'Measured');

% Plot model curve as dashed line
plot(freq_plot, y_fit, 'r--', 'LineWidth', 2, 'DisplayName', 'Fitted Model');


xlabel('Frequency [Hz]', 'FontSize', 16);
ylabel('Amplitude ratio |V_R/V_{in}|', 'FontSize', 16);
title(sprintf('Amplitude Ratio Fit (\\chi^2_{red} = %.3f)', reduced_chi2), 'FontSize', 16);

legend('Location', 'best');
set(gca, 'FontSize', 14);

% Print parameter values and their errors to the console
fprintf('Fitted Parameters:\n');
fprintf('L/R   = %.3e ± %.1e\n', params_fit(1), param_errors(1));
fprintf('1/RC  = %.3e ± %.1e\n', params_fit(2), param_errors(2));



% Convert frequencies to angular frequencies
w_R = 2 * pi * freq_H(:);  % Angular frequency vector [rad/s]

% Model: phi(w) = atan(L/R * w - 1/(RC * w))
modelPhase = @(p, w) atan(-p(1) * w + p(2) ./ w);  % p(1) = L/R, p(2) = 1/RC

% Initial guess based on expected component values
startPoint = [H / 470, 1 / (470 * C)];  % [L/R, 1/RC]

% Bounds
lb = [0, 0];
ub = [Inf, Inf];

% Target data
yData = phase_diff(:)*-1-pi/2;

% Perform nonlinear least squares fitting
[p_fit, resnorm, residuals, exitflag, output, lambda, J] = lsqcurvefit( ...
    modelPhase, startPoint, w_R, yData, lb, ub);

% Confidence interval estimation
ci = nlparci(p_fit, residuals, 'jacobian', J);
param_errors = (ci(:,2) - ci(:,1)) / 2;

% Compute reduced chi-squared
N = length(yData);      % number of data points
k = length(p_fit);      % number of fitted parameters
y_fit_data = modelPhase(p_fit, w_R);
chi2_red = sum(((yData - y_fit_data) ./ err_phase(:)).^2) / (N - k);

% Display results
fprintf('Phase Fit Parameters:\n');
fprintf('L/R   = %.3e ± %.1e\n', p_fit(1), param_errors(1));
fprintf('1/RC  = %.3e ± %.1e\n', p_fit(2), param_errors(2));
fprintf('Reduced Chi^2 = %.3f\n', chi2_red);

% Plot results
figure; hold on;
errorbar(freq_H, yData, err_phase(:), '.', 'MarkerSize', 25, 'DisplayName', 'Measured');
w_plot = linspace(min(w_R), max(w_R), 500)';
f_plot = w_plot / (2*pi);  % Convert back to Hz for x-axis
y_fit = modelPhase(p_fit, w_plot);
plot(f_plot, y_fit, 'r--', 'LineWidth', 2, 'DisplayName', 'Fitted Model');

xlabel('Frequency [Hz]', 'FontSize', 16);
ylabel('Phase diff [rad]', 'FontSize', 16);
title('phase difference as a function of input frequency')
legend('Location', 'best');
set(gca, 'FontSize', 14);

clear alpha beta cov_matrix delta dx Fs gamma i jacobian k L lb model N P paramsFitted1 paramsFitted2 refinedFreq reifnedIndex residual std_errors_1 std_errors_2 ub var_res w xdata Y y_detrended y_windowed ydata B0 data1 data2
