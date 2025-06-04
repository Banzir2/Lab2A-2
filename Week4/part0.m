close all; clear;

data1 = readmatrix('data/work_of_art_squared/meas1.csv');
data2 = readmatrix('data/work_of_art_squared/meas2.csv');

figure; hold on;
scatter(data1(:, 1), data1(:, 2));
scatter(data2(:, 1), data2(:, 2));

