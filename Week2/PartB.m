close all; clear;
C_33n=37.45*1e-9;
L=2.152/1000;
R_100k=100.1*1e3;
R_1k=999;
freq1k = [1.143*1e3, 111.15, 120.7*1e3, 128.4, 157.2, 200.4,...
    22.77, 239.9, 25.11*1e3, 44.70, 5.052*1e3, 68.03, 89.61];
freq100k = [1.309*1e3, 107.3, 120.2*1e3, 131.4, 165.3, 191.6,...
    223.2, 24.33*1e3, 244.9, 37.16, 5.340*1e3, 59.31, 79.11];

file_names = ["1.1k-1k", "1.3k-100k", "107-100k", "111-1k", "120k-100k", "121k-1k", "128-1k", "131-100k", ...
                      "157-1k", "165-100k", "191-100k", "200-1k", "22-1k", "223-100k", "239-1k", "24.3k-100k", ...
                      "245-100k", "25k-1k", "37-100k", "44-1k", "5.3k-100k", "59-100k", "5k-1k", "68-1k", "79-100k", ...
                      "89-1k"];
indexes100k = [2, 3, 5, 8, 10, 11, 14, 16, 17, 19, 21, 22, 25];
indexes1k = [1, 4, 6, 7, 9, 12, 13, 15, 18, 20, 23, 24, 26];

ampl_diff = zeros(size(file_names));
phase_diff = zeros(size(file_names));
err_amp = zeros(size(file_names));
err_phase = zeros(size(file_names));
for i = 1:size(file_names, 2)
    data1 = readmatrix(strcat('PartB/', file_names(i), '/meas1.csv'));
    data2 = readmatrix(strcat('PartB/', file_names(i), '/meas2.csv'));

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

amp_diff_1k = ampl_diff(indexes1k);
amp_diff_100k = ampl_diff(indexes100k);
phase_diff_1k = phase_diff(indexes1k);
phase_diff_100k = phase_diff(indexes100k);
err_amp_1k = err_amp(indexes1k);
err_amp_100k = err_amp(indexes100k);
err_phase_1k = err_phase(indexes1k);
err_phase_100k = err_phase(indexes100k);

xhisquare = zeros(1, 4);

figure; hold on; % Amp diff 1k ohm
errorbar(log(freq1k), amp_diff_1k, err_amp_1k * 4, '.', 'MarkerSize', 25);
f = fit(freq1k.', amp_diff_1k.', '(x / sqrt(x^2 + a)) + b');
diff2 = (amp_diff_1k - f(freq1k).').^2;
xhisquare(1) = sum(diff2 ./ (err_amp_1k.^2)) / 10;
xlabel('Log(Frequency)', 'FontSize', 16);
ylabel('Amplitude ratio', 'FontSize', 16);
title('Amplitude ratio by input frequency 1K\Omega', 'FontSize', 16);
ax = gca;
ax.FontSize = 14;

figure; hold on; % Phase diff 1k ohm
errorbar(log(freq1k), phase_diff_1k, err_phase_1k * 30, '.', 'MarkerSize', 25);
f = fit(freq1k.', phase_diff_1k.', 'atan(a/x)');
diff2 = (phase_diff_1k - f(freq1k).').^2;
xhisquare(2) = sum(diff2 ./ (err_phase_1k.^2 + (pi*freq1k*1e-4).^2)) / 10;
xlabel('Log(Frequency)', 'FontSize', 16);
ylabel('Phase diff [rad]', 'FontSize', 16);
title('Phase diff by input frequency 1K\Omega', 'FontSize', 16);
ax = gca;
ax.FontSize = 14;

figure; hold on; % Amp diff 100k ohm
errorbar(log(freq100k), amp_diff_100k, err_amp_100k, '.', 'MarkerSize', 25);
f = fit(freq100k.', amp_diff_100k.', 'x / sqrt(x^2 + a^2)');
diff2 = (amp_diff_100k - f(freq100k).').^2;
xhisquare(3) = sum(diff2 ./ (err_amp_100k.^2 )) / 10;
xlabel('Log(Frequency)', 'FontSize', 16);
ylabel('Amplitude ratio', 'FontSize', 16);
title('Amplitude ratio by input frequency 100K\Omega', 'FontSize', 16);
ax = gca;
ax.FontSize = 14;

figure; hold on; % Phase diff 100k ohm
errorbar(log(freq100k), phase_diff_100k, err_phase_100k * 30, '.', 'MarkerSize', 25);
f = fit(freq100k.', phase_diff_100k.', 'atan(a/x) + b');
diff2 = (phase_diff_100k - f(freq100k).').^2;
xhisquare(4) = sum(diff2 ./ (err_phase_100k.^2 + (pi*freq1k*1e-4).^2)) / 10;
xlabel('Log(Frequency)', 'FontSize', 16);
ylabel('Phase diff [rad]', 'FontSize', 16);
title('Phase diff by input frequency 100K\Omega', 'FontSize', 16);
ax = gca;
ax.FontSize = 14;