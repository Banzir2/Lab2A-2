close all; clear;
freq_H=[33.8,36.39, 38.58, 41.25, 45.17, 31.13, 29.38, 27.29, 23.19, 105.7, 10.15]*1e3;
file_names = ["34k", "36k", "38k", "41k", "45k", "31k", "27k", "29k", "23k","106k","10k"];
xhisquare = zeros(1, 2);
C = 1.0571e-09;
H = 21.06e-3;
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
p_fit_H = p_fit;
param_errors_H = param_errors;
w_H = w_R;
modelPhase_H = modelPhase;
phase_diff_H = phase_diff;
err_phase_H = err_phase;

% Frequencies in Hz
freq_R = [35.66, 37.54, 39.37, 46.99, 33.03, 31.25, 28.12, 24.13, 100.7, 10.02]*1e3;
file_names = ["35k", "37k", "39k", "47k", "33k", "31k", "28k", "24k", "100k","10k"];
w = 2 * pi * freq_R(:);    % Angular frequencies

% Initialize
ampl_diff = zeros(size(file_names));
phase_diff = zeros(size(file_names));
err_amp = zeros(size(file_names));
err_phase = zeros(size(file_names));

for i = 1:length(file_names)
    data1 = readmatrix(strcat('470ohm-b/', file_names(i), '/meas1.csv'));
    data2 = readmatrix(strcat('470ohm-b/', file_names(i), '/meas2.csv'));

    % Model function to fit: y = a * exp(b * x)
    freq = freq_R(i);
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
% Remove bad data point(s)
bad_idx = 2;  % Adjust if needed
phase_diff(bad_idx) = [];
freq_R(bad_idx) = [];
err_phase(bad_idx) = [];

% Convert frequencies to angular frequencies
w_R = 2 * pi * freq_R(:);  % Angular frequency vector [rad/s]

% Model: phi(w) = atan(L/R * w - 1/(RC * w))
modelPhase = @(p, w) atan(-p(1) * w + p(2) ./ w);  % p(1) = L/R, p(2) = 1/RC

% Initial guess based on expected component values
startPoint = [25e-3 / 470, 1 / (470 * 1e-9)];  % [L/R, 1/RC]

% Bounds
lb = [0, 0];
ub = [Inf, Inf];

% Target data
yData = phase_diff(:)*(-1);

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
errorbar(freq_R, yData, err_phase(:), '.', 'MarkerSize', 25, 'DisplayName', 'Measured');
w_plot = linspace(min(w_R), max(w_R), 500)';
f_plot = w_plot / (2*pi);  % Convert back to Hz for x-axis
y_fit = modelPhase(p_fit, w_plot);
plot(f_plot, y_fit, 'r--', 'LineWidth', 2, 'DisplayName', 'Fitted Model');

xlabel('Frequency [Hz]', 'FontSize', 16);
ylabel('Phase diff [rad]', 'FontSize', 16);
title('phase difference as a function of input frequency')
legend('Location', 'best');
set(gca, 'FontSize', 14);

p_fit_R = p_fit;
param_errors_R = param_errors;
modelPhase_R = modelPhase;
phase_diff_R = phase_diff;
err_phase_R = err_phase;

R=[10,100,220,470, 1000];
err_R=[13.8,99.6,220.3,465,1004];

freq_10=[32.98, 36.5, 30.05];
Amplitude_10ohm=[18.6,6.24,6.24];
Amplitude_otot10ohm=[1.1, 1.26, 1.28];

freq_100=[32.97,37.23,29.08];
Amplitude_otot100ohm=[1.14,1.26, 1.29 ];
Amplitude_100ohm=[14.8, 5.04, 4.88];

freq_1k=[33.2,44.49, 13.93];
Amplitude_otot1kohm=[1.26,1.3, 1.3];
Amplitude_1kohm=[4.56, 1.52,1.52];

freq_220=[33.07,38.45, 27.74];
Amplitude_otot10kohm=[1.18,1.28, 1.29];
Amplitude_10kohm=[11.6, 3.8,3.84];

freq_420=[35.11, 37.37, 39.25, 43.18, 46.38, 33.16, 31.29, 28.03, 24.34, 20.18, 83.19, 168.4, 10.86, 5.215]*1e3;

file_names = ["35k", "37k", "39k", "43k", "46k", "33k", "31k", "28k", "24k", "20k", "83k", "168k", "10k", "5k"];

ampl_diff = zeros(size(file_names));
phase_diff = zeros(size(file_names));
err_amp = zeros(size(file_names));
err_phase = zeros(size(file_names));
for i = 1:size(file_names, 2)
    data1 = readmatrix(strcat('470ohm/', file_names(i), '/meas1.csv'));
    data2 = readmatrix(strcat('470ohm/', file_names(i), '/meas2.csv'));

    % Model function to fit: y = a * exp(b * x)
    freq = freq_420(i);
    model = @(params, x) params(1) * cos(2*pi*freq * x + params(2));

    lb = [0, -pi];
    ub = [100, pi];

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

[freq_420, newIdx] = sort(freq_420);       % I gives the old indices in new order
ampl_diff = ampl_diff(newIdx);
phase_diff = phase_diff(newIdx);
err_amp = err_amp(newIdx);
err_phase = err_phase(newIdx);

xhisquare = zeros(1, 2);
C = 1.0571e-09;
H = 21.06e-3;
freq_space = linspace(0, max(freq_420)*1.01, 10000);
figure; hold on; % Phase diff
errorbar(freq_420, -(phase_diff - pi/2), err_phase * 30, '.', 'MarkerSize', 25);
modelFun = @(params, x) atan((params(1) ./ x) - (x / params(2)));
startPoint = [1 / (2*pi*470*C), 470 / (2*pi*H)];
xData = freq_420(:);     % Ensure column vector
yData = -(phase_diff.' - pi/2);
lb = [0, 0];
ub = [Inf, Inf];
options = optimoptions('lsqcurvefit', 'TolFun', 1e-8);
[params_fit, ~, residual, ~, ~, ~, jacobian] = lsqcurvefit(modelFun, startPoint, xData, yData, lb, ub, options);
var_res = sum(residual.^2) / (length(residual) - length(params_fit));
cov_matrix = var_res * inv(jacobian' * jacobian);
std_errors = full(sqrt(diag(cov_matrix)));
text(0.7 * max(freq_space), 0.5*max(phase_diff), compose("f_c = %.1f \\pm %.1f", params_fit(1), std_errors(1)), 'FontSize', 16, 'FontName', 'Times New Roman');
text(0.7 * max(freq_space), 0.5*max(phase_diff) - 0.25, compose("f_l = %.1f \\pm %.1f", params_fit(2), std_errors(2)), 'FontSize', 16, 'FontName', 'Times New Roman');
plot(freq_space, modelFun(params_fit, freq_space), 'LineWidth', 1.5);
diff2 = (phase_diff + modelFun(params_fit, freq_420) - pi/2).^2;
xhisquare(2) = sum(diff2 ./ (err_phase.^2 + (pi*freq_420*1e-7).^2)) / 12;
xlabel('Frequency [Hz]', 'FontSize', 16);
ylabel('Phase diff [rad]', 'FontSize', 16);
title('Phase diff by input frequency', 'FontSize', 16);
ax = gca;
ax.FontSize = 14;

% === MERGED PLOT OF ALL THREE PHASE FITS ===
figure; hold on;

% Plot 1: High-capacitance set (H)
errorbar(freq_H, -phase_diff_H - pi/2, err_phase_H, '.', 'MarkerSize', 18, ...
    'HandleVisibility', 'off', 'Color', [1, 0.4, 0.4]);

% Plot fitted model (H)
w_plot_H = linspace(min(w_H), max(w_H), 500)';
f_plot_H = w_plot_H / (2*pi);
y_fit_H = modelPhase_H(p_fit_H, w_plot_H);
plot(f_plot_H, y_fit_H, '--', 'LineWidth', 2, 'Color', [1, 0.4, 0.4], ...
    'DisplayName', 'inductor');

% Plot 2: R capacitor
errorbar(freq_R, -phase_diff_R, err_phase_R, '.', 'MarkerSize', 18, ...
    'HandleVisibility', 'off', 'Color', [0.1, 0.7, 0.1]);

w_plot_R = linspace(min(w_R), max(w_R), 500)';
f_plot_R = w_plot_R / (2*pi);
y_fit_R = modelPhase_R(p_fit_R, w_plot_R);
plot(f_plot_R, y_fit_R, '--', 'LineWidth', 2, 'Color', [0.1, 0.7, 0.1], ...
    'DisplayName', 'resistor');

% Plot 3: General 470 Ohm
errorbar(freq_420, -(phase_diff - pi/2), err_phase , '.', 'MarkerSize', 18, ...
    'HandleVisibility', 'off', 'Color', [0.2, 0.6, 1]);

y_fit_general = modelFun(params_fit, freq_space);
plot(freq_space, y_fit_general, '--', 'LineWidth', 2, 'Color', [0.2, 0.6, 1], ...
    'DisplayName', 'capicator');

% === Formatting ===
xlabel('Frequency [Hz]', 'FontSize', 16);
ylabel('Phase difference [rad]', 'FontSize', 16);
title('Phase difference vs Input Frequency', 'FontSize', 18);
legend('Location', 'best');
set(gca, 'FontSize', 14);
grid on;
