% in order to calculate the moment arm we extracted data from OpenSim. We
% use the data to find the moment arm at 90 degree
clear all
close all
clc
% import the data
importdata('momentarm.txt');

m1 = ans.data(:,2);
x  = ans.data(:,3);

% averaging two points (one before and one after the desired point)
xmarkers = (90.9091+89.5455)/2;
ymarkers = (0.0517+0.0522)/2;

% plot the data
plot(m1,x, 'linewidth', 2)
hold on
plot(xmarkers, ymarkers, 'r*');set(gca,'fontsize',14);
xlabel('Elbow Angle [degree]');
ylabel('Moment Arm [m]');

% define the desired point
x1 = xmarkers;
y1 = ymarkers;
str1 = '    M at 90 = 0.052';
text(x1, y1, str1);set(gca,'fontsize',14);
