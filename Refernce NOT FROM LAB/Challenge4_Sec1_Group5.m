%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CODE CHALLENGE 4 - Linear Least-Squares Fit
%
% The purpose of this program is to calculate the equation of the best fit
% line for a data set using linear least-squares fitting.
%
% To complete the challenge, finish the code below to:
% 1) load data from csv file
% 2) find linear best fit coefficients and associated uncertainty
% 3) plot the original data along with the best fit line 
% 4) add errorbars for fit uncertainty to this plot from the data and from
%    the linear regression parameters
%
% NOTE: DO NOT change any variable names already present in the code.
% 
% Upload your team's script to Gradescope when complete.
% 
% NAME YOUR FILE AS Challenge4_Sec{section number}_Group{group breakout #}.m 
% ***Section numbers are 1 or 2*** 
% EX File Name: Challenge4_Sec1_Group15.m 
%
% STUDENT TEAMMATES (Group 5)
% 1) AJ Lauffer
% 2) Skylar Clark skcl3151@colorado.edu
% 3) Hattie Rice hari2695@colorado.edu
% 4) Muhannad W Ibrahim muib1439@colorado.edu
% 5) Ceu Gomez-Faulk cego6160@colorado.edu

% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Housekeeping (Please don't "clear all" or "clearvars", it makes grading difficult)
close all   % Close all open figure windows
clc           % Clear the command window

%% Load and extract the time and velocity vectors from the data
data = xlsread("Challenge4_data.csv");
t = data(:,1);    % [s]
v = data(:,2);    % [m/s]

%% Calculations
% Find number of data points in the vectors
N = length(t);

% Find linear best fit coefficients A and B
% Create H matrix
H = [ones(length(t),1) t];

% Create y matrix
y = v;

% Create W matrix (hint: type <help diag> in command line)
W = diag(ones(1,N)); % set W as identity matrix
SigmaY = std(v)/sqrt(N-2); % compute SigmaY from Eqn 8.15
W = (1/SigmaY)*W;% recompute the W matrix using SigmaY


% Solve for P matrix
P = (H' * W * H)^(-1);

% Solve for x_hat matrix and extract A and B parameters
x_hat = P * H' * W * y;
A = x_hat(1);
B = x_hat(2);

% extract uncertainty in A and uncertainty in B from P matrix
A_error = sqrt(P(1));
B_error = sqrt(P(4));

%% Display acceleration with associated uncertainty and the intial velocity with associated uncertainty
%  Make sure to use and display with CORRECT SIG FIGS
fprintf("Acceleration is %.1f ± %.1f \n",B,B_error);
fprintf("Inital velocity is %.1f ± %.1f \n",A,A_error);

%% Find predicted velocity values using your linear fit equation
v_predicted = A + B * t;

%% Ploting and Error Calculations
% On the same plot, do the following:
% 1. plot the velocity data vs time as a scatter plot 
hold on
scatter(t,v);
% 2. plot predicted velocity vs time as a line
plot(t,v_predicted);
% 3. title your plot so that it indicates what celestial body this data
%    simulates
title("Free Fall on the Moon");
% 4. Add measured velocity error bars and predicted velocity error bars to 
%    the plot (hint - this will involve error propagation calculations
v_err = std(v)/sqrt(N);
v_predicted_error = H * P * H';
errorbar(t,v,v_err*ones(N,1));
errorbar(t,v,v_predicted_error*ones(N,1));
xlabel("Time (s)");
ylabel("Velocity (m/s)");
hold off