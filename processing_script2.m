% Read the CSV file
clc, clearvars
clear all

disp(pwd);
addpath('C:/Users/ihsan/OneDrive/Documents/coding/ECE150/Other_CPP/Lorentz-Solver');
data = readtable("data.csv");

% Extract the data
x = data.x;
y = data.y;
z = data.z;

% Create a 3D plot
figure;
plot3(x, y, z, '-o');

% Add labels and title
xlabel('X');
ylabel('Y');
zlabel('Z');
title('3D Plot of CSV Data');
grid on;