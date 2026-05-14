clc;
clear;
close all
%%
f_c = 3e9;
SNR = 20;
interval =8; % pilot interval
N=256;  % number of subcarriers
M=64;  %number of antennas
G=20 ; % number of frames
L = 30; % number of delay taps
bandwidth = 2e7;
%%
[Y, H, sequence_Phi, coefficients,F_truncated] = signal_generation(SNR, interval, N, G, M, L, bandwidth);
H_est = EPA(Y, M,  coefficients,sequence_Phi,F_truncated);
norm( H - H_est, 'fro' )^2 / norm( H, 'fro' )^2



