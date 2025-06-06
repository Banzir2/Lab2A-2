close all; clear;
R=[10,100,220,470, 1000];
err_R=[13.8,99.6,220.3,465,1004];

freq_10=[32.98, 36.5, 30.05]*1e3;
Amplitude_10ohm=[18.6,6.24,6.24];
Amplitude_otot10ohm=[1.1, 1.26, 1.28];

freq_100=[32.97,37.23,29.08]*1e3;
Amplitude_otot100ohm=[1.14,1.26, 1.29 ];
Amplitude_100ohm=[14.8, 5.04, 4.88];

freq_1k=[33.2,44.49, 13.93]*1e3;
Amplitude_otot1kohm=[1.26,1.3, 1.3];
Amplitude_1kohm=[4.56, 1.52,1.52];

freq_220=[33.07,38.45, 27.74]*1e3;
Amplitude_otot220ohm=[1.18,1.28, 1.29];
Amplitude_220ohm=[11.6, 3.8,3.84];

freq_420=[35.11, 37.37, 39.25, 43.18, 46.38, 33.16, 31.29, 28.03, 24.34, 20.18, 83.19, 168.4, 10.86, 5.215]*1e3;

file_names = ["35k", "37k", "39k", "43k", "46k", "33k", "31k", "28k", "24k", "20k", "83k", "168k", "10k", "5k"];

ampl_diff = zeros(size(file_names));
C_amp = zeros(size(file_names));
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

    C_amp(i) = abs(paramsFitted2(1));
    ampl_diff(i) = abs(paramsFitted2(1) / paramsFitted1(1));
    phase_diff(i) = mod(paramsFitted1(2) - paramsFitted2(2), pi) - pi;
    err_amp(i) = scope_err(data1(1, 3)) + scope_err(data2(1, 3));
    err_phase(i) = std_errors_1(2) + std_errors_2(2);
end

clear alpha beta cov_matrix delta dx Fs gamma i jacobian k L lb model N P paramsFitted1 paramsFitted2 refinedFreq reifnedIndex residual std_errors_1 std_errors_2 ub var_res w xdata Y y_detrended y_windowed ydata B0 data1 data2

new_R = ones(1, 6) * 470;
new_R_err = ones(1, 6);

[freq_420, newIdx] = sort(freq_420);       % I gives the old indices in new order
ampl_diff = ampl_diff(newIdx);
phase_diff = phase_diff(newIdx);
err_amp = err_amp(newIdx);
err_phase = err_phase(newIdx);

xhisquare = zeros(1, 2);
C = 1.0571e-09;
H = 21.06e-3;

figure; hold on; % Amp diff
errorbar(freq_420, ampl_diff, err_amp, '.', 'MarkerSize', 25);

modelFun = @(params, x) params(1) ./ (sqrt(1 + (x/params(2) - params(1)./x).^2) .* x);
startPoint = [1 / (2*pi*470*C), 470 / (2*pi*H)];
xData = freq_420(:);     % Ensure column vector
yData = ampl_diff(:);
lb = [0, 0];
ub = [Inf, Inf];
options = optimoptions('lsqcurvefit', 'TolFun', 1e-8);
[params_fit, ~, residual, ~, ~, ~, jacobian] = lsqcurvefit(modelFun, startPoint, xData, yData, lb, ub, options);
var_res = sum(residual.^2) / (length(residual) - length(params_fit));
cov_matrix = var_res * inv(jacobian' * jacobian);
std_errors = full(sqrt(diag(cov_matrix)));

new_R(1) = 1 / (2*pi*C*params_fit(1)) - 470;
new_R_err(1) = (1 / (2*pi*C*(params_fit(1) - std_errors(1))) - 1 / (2*pi*C*(params_fit(1) + std_errors(1)))) / 2;
new_R(2) = 2*pi*H*params_fit(2) - 470;
new_R_err(2) = 2*pi*H*std_errors(2);
freq_space = linspace(0, max(freq_420)*1.01, 10000);

plot(freq_space, modelFun(params_fit, freq_space), 'LineWidth', 1.5);
text(0.7 * max(freq_space), max(ampl_diff), compose("f_c = %.1f \\pm %.1f", params_fit(1), std_errors(1)), 'FontSize', 16, 'FontName', 'Times New Roman');
text(0.7 * max(freq_space), max(ampl_diff) - 0.4, compose("f_l = %.1f \\pm %.1f", params_fit(2), std_errors(2)), 'FontSize', 16, 'FontName', 'Times New Roman');
diff2 = (ampl_diff - modelFun(params_fit, freq_420)).^2;
xhisquare(1) = sum(diff2 ./ (err_amp.^2)) / 12;
xlabel('Frequency [Hz]', 'FontSize', 16);
ylabel('Amplitude ratio', 'FontSize', 16);
title('Amplitude ratio by input frequency', 'FontSize', 16);
xline(4.004209208122798e+04, 'LineWidth', 2, 'LineStyle', '--');
xline(2.492878536140029e+04, 'LineWidth', 2, 'LineStyle', '--');
ax = gca;
ax.FontSize = 14;

thirdmax = max(modelFun(params_fit, freq_space)) / 3;
fun = @(x) modelFun(params_fit, x) - thirdmax;
x = fzero(fun, 6e4);

figure; hold on; % Phase diff
errorbar(freq_420, phase_diff, err_phase * 30, '.', 'MarkerSize', 25);
modelFun = @(params, x) atan((x / params(2)) - (params(1) ./ x)) - pi/2;
startPoint = [1 / (2*pi*470*C), 470 / (2*pi*H)];
xData = freq_420(:);     % Ensure column vector
yData = phase_diff.';
lb = [0, 0];
ub = [Inf, Inf];
options = optimoptions('lsqcurvefit', 'TolFun', 1e-8);
[params_fit, ~, residual, ~, ~, ~, jacobian] = lsqcurvefit(modelFun, startPoint, xData, yData, lb, ub, options);
var_res = sum(residual.^2) / (length(residual) - length(params_fit));
cov_matrix = var_res * inv(jacobian' * jacobian);
std_errors = full(sqrt(diag(cov_matrix)));

new_R(3) = 1 / (2*pi*C*params_fit(1)) - 470;
new_R_err(3) = (1 / (2*pi*C*(params_fit(1) - std_errors(1))) - 1 / (2*pi*C*(params_fit(1) + std_errors(1)))) / 2;
new_R(4) = 2*pi*H*params_fit(2) - 470;
new_R_err(4) = 2*pi*H*std_errors(2);

text(0.7 * max(freq_space), max(phase_diff) - 0.5, compose("f_c = %.1f \\pm %.1f", params_fit(1), std_errors(1)), 'FontSize', 16, 'FontName', 'Times New Roman');
text(0.7 * max(freq_space), max(phase_diff) - 0.75, compose("f_l = %.1f \\pm %.1f", params_fit(2), std_errors(2)), 'FontSize', 16, 'FontName', 'Times New Roman');
plot(freq_space, modelFun(params_fit, freq_space), 'LineWidth', 1.5);
diff2 = (phase_diff - modelFun(params_fit, freq_420)).^2;
xhisquare(2) = sum(diff2 ./ (err_phase.^2 + (pi*freq_420*1e-7).^2)) / 12;
xlabel('Frequency [Hz]', 'FontSize', 16);
ylabel('Phase diff [rad]', 'FontSize', 16);
title('Phase diff by input frequency', 'FontSize', 16);
ax = gca;
ax.FontSize = 14;

f_c = @(r) 1 ./ (2*pi*r*C);
f_l = @(r) r ./ (2*pi*H);
a = @(r) 1 ./ (f_l(r).^2);
b = @(r) 1 - 2*(f_c(r) ./ f_l(r));
c = @(r) f_c(r).^2 - 9*f_c(r).*f_l(r);
w = @(param, r) sqrt((-b(r+param(1)) + sqrt(b(r+param(1)).^2 - 4*a(r+param(1)).*c(r+param(1)))) ./ (2*a(r+param(1)))) - sqrt((-b(r+param(1)) - sqrt(b(r+param(1)).^2 - 4*a(r+param(1)).*c(r+param(1)))) ./ (2*a(r+param(1))));
 
figure; hold on;
R = [10, 100, 220, 470, 1000];
widths = [freq_10(2) - freq_10(3), freq_100(2) - freq_100(3), ...
    freq_220(2) - freq_220(3), 1.511330671982769e+04, freq_1k(2) - freq_1k(3)];
errorbar(R, widths, 5*widths*1e-2, '.', 'MarkerSize', 25);
startPoint = [2.538035872248687e+02];
lb = [0];
ub = [1000];
options = optimoptions('lsqcurvefit', 'TolFun', 1e-8);
[params_fit, ~, residual, ~, ~, ~, jacobian] = lsqcurvefit(w, startPoint, R', widths', lb, ub, options);
var_res = sum(residual.^2) / (length(residual) - length(params_fit));
cov_matrix = var_res * inv(jacobian' * jacobian);
std_errors = full(sqrt(diag(cov_matrix)));

new_R(5) = params_fit;
new_R_err(5) = std_errors;

xhi = sum((residual' ./ (5*widths*1e-2)).^2) / 4;

plot(0:0.1:1000, w(params_fit, 0:0.1:1000), 'LineWidth', 1.5);
xlabel('Resistance [ohm]', 'FontSize', 16);
ylabel('Curve Width [Hz]', 'FontSize', 16);
title('Curve Width by Resistance', 'FontSize', 16);
ax = gca;
ax.FontSize = 14;

figure; hold on;
amps = [Amplitude_10ohm(1) / Amplitude_otot10ohm(1), Amplitude_100ohm(1) / Amplitude_otot100ohm(1), Amplitude_220ohm(1) / Amplitude_otot220ohm(1), max(ampl_diff), Amplitude_1kohm(1) / Amplitude_otot1kohm(1)];
scatter(R, amps, '.', 'SizeData', 600);
f = fit(R.', amps.', '4.4635e+03/(x+a)', 'StartPoint', 220);
new_R(6) = f.a;
new_R_err(6) = 10;
plot(f);

to_print = ["a", "a", "a", "a", "a", "a"];
for i=1:6
    fprintf(compose("%f+-%f\n", new_R(i), new_R_err(i)));
end