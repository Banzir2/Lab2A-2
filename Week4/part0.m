close all; clear;

data1 = readmatrix('data/work_of_art_squared/meas1.csv');
data2 = readmatrix('data/work_of_art_squared/meas2.csv');

figure; hold on;
plot(data1(:, 1), data1(:, 2), 'LineWidth', 2);
plot(data2(:, 1), fliplr(data2(:, 2)')', 'LineWidth', 2);
% plot(data1(:, 1), data1(:, 2) + data2(:, 2), 'LineWidth', 2);
xlabel('Time [sec]', 'FontSize', 16);
ylabel('Voltage [volt]', 'FontSize', 16);
title('R1 and R2 voltage by time (8V input)', 'FontSize', 16);
legend('R1', 'R2', 'R1 + R2', 'FontSize', 16);
ax = gca;
ax.FontSize = 14;