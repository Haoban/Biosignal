function y = chirpmodel(t, delay, a, b, lambda_w, f0, f1, lambda_f, model)
% Copyright (c) 2015, Eero Väyrynen.
% CHIRPMODEL.m generates linear, 3rd order and 4th order polynomial approximation,
% exponential, and linear to exponential amplitude dampened single instantaneous 
% frequency chirps using the supplied parameters (and global variables)
% with an initial phase of pi/2 (cos -> sin). 
%
%
% INPUTS
%
% t                         Linear time vector of desired length
% delay                     Onset delay of the chirp
% a                         Initial amplitude weight
% b                         End amplitude weight
% lambda_w                  Amplitude nonlinearity coefficient
% f0                        Initial frequency
% f1                        Target end frequency
% lambda_f                  Frequency nonlinearity coefficient
% model                     Chirp model (lin, poly3, poly4, exp, linexp)
%
% OUTPUTS
%
% y                         Chirp waveform

global fs                   % Sampling frequency
global chirp_length         % Chirp length in samples
delayed_t = t-delay;        % Apply delay
dur = chirp_length/fs;      % Chirp duration in seconds

switch model
    case 'lin'
    chirp_y = sin(2*pi*(f0.*delayed_t + ((f1-f0)/(dur*2)).*delayed_t.^2)); % linear chirp

    case 'poly3'
    chirp_y = sin(2*pi*(f0.*delayed_t + ((f1-f0)/(dur*2)).*delayed_t.^2 + (f1-f0)*(1-lambda_f)*(1/(dur*4).*delayed_t.^2 - 1/(dur^(2)*6).*delayed_t.^(3)))); % 3rd order polynomial approximation chirp
        
    case 'poly4'
    chirp_y = sin(2*pi*(f0.*delayed_t + ((f1-f0)/(dur*2)).*delayed_t.^2 + (f1-f0)*(1-lambda_f)*(2/(dur*2).*delayed_t.^2 - 6*(1/(dur^(2)*6)).*delayed_t.^(3) + 6*(1/(dur^(3)*24)).*delayed_t.^(4)))); % 4th order polynomial approximation chirp
        
    case 'exp'
    k = f1/f0;
    chirp_y = sin(((dur*2*pi*f0)/log(k)).*(exp(delayed_t./(dur/log(k)))-1)); % exponential chirp

    case 'linexp'
    k = f1/f0;
    chirp_y = sin((1-lambda_f)*((dur*2*pi*f0).*(k.^(delayed_t/dur)-1)/log(k))+(lambda_f)*(2*pi*f0.*delayed_t + 2*pi*((f1-f0)/(dur*2)).*delayed_t.^2)); % exp + linear

    otherwise
        error('No model or unknown model specified!')
end

A = (a-b).*(exp(-delayed_t*lambda_w)-exp(-dur*lambda_w))./(1-exp(-dur*lambda_w))+b; % Linear to exponential dampening function

y = chirp_y.*A; % apply dampening
y( delayed_t < 0 ) = 0; % set signal to zero becore chirp onset


