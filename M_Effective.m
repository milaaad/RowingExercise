% Calculating effective mass of the forearm and hand which is used in the
% equation of motion as m3
% All the data are extracted from Winter's table
% To calculate the effective mass we need to know mass, location of the
% center of mass, and moment of inertia
% Winter's data table may be found online (http://www.ele.uri.edu/faculty/vetter/BME207/anthropometric-data.pdf)
clear all
clc
Mbody = 75;                 % (kg) body mass
Marm  = 0.022*Mbody;        % (kg) Forearm + hand mass
Larm  = 0.517;              % (m) Length of the forearm+hand 
RG    = 0.827*Larm;         % (m) Radius of Gyration for forearm+hand, with respect to elbow
Iarm  =  Marm*RG^2;         % (kg*m^2) Moment of inertia

Meff  = Iarm / (Larm^2);
fprintf('Effective Mass = %8.4f kg			\n', Meff);