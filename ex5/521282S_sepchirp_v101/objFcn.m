function [out] = objFcn(in)
% Copyright (c) 2015, Eero Väyrynen, version 1.01
% objFUN.m is an objective fuction used by the Particle Swarm Optimization. The
% function implements a linear weight of errors to emphasise the beginnings
% of the chirp models.
%
% Inputs (in) are N particle length vectors of parameter values. Other
% parameters used by the objective function are stored in the global variables.
%
% Outputs (out) are single objective function values for each particle in
% vector of length N

% Version 1.01 changes
% - linear MSE calculation bug fix (corrected erroneous nonlinear weight)

global fs               % sampling frequency of the data
global data_points      % length of the data sample
global origData         % Original data used to test the chirp models
global chirp_length     % length of the chirp generated
global model            % chirp model

delay = in(:,1);        % Chirp onset delay (Tau)
a = in(:,2);            % Amplitude initial weight (a)
b = in(:,3);            % Amplitude end weigth (b)
lambda_w = in(:,4);     % Amplitude nonlinearity (lambda_w)
f0 = in(:,5);           % Start frequency (f0)
f1 = in(:,6);           % End frequency (f1)
lambda_f = in(:,7);     % Frequency nonlinearity (lambda_f)

G = zeros(length(a),1); % preallocate G
chirp_t = linspace(0,data_points/fs,data_points); % Linear time vector
N = length(a); % get number of particles

for i = 1:N,
    optChirp = chirpmodel(chirp_t, delay(i), a(i), b(i), lambda_w(i), f0(i), f1(i), lambda_f(i), model); % call chirpmodel

    err = (optChirp-origData); %error
    
    % Linear emphasis of the chirp part
    delay_zeros = floor(delay(i)*fs); % number of initial zeros
    se = (err.^2).*[linspace(1,1,delay_zeros) linspace(1,0.1,chirp_length) linspace(0.1,0.1,(data_points-delay_zeros-chirp_length))]; % squared error with linear weight
    
    % Mean squared error 
    mse = mean(se);

    G(i) = mse; % collect mse to vector output
end

out = G;

end
