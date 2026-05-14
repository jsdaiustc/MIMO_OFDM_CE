function [Y, H_true, sequence_Phi, coefficients ,F_truncated] = signal_generation(SNR, pilot_interval, N, G, M, L, bandwidth)
pilot_inv = pilot_interval;
B = bandwidth;
pilot_loc = [];
data_loc = [];
for k = 1:N
    if(rem(k,pilot_inv)~=1)
        data_loc = [data_loc k];
    else
        pilot_loc = [pilot_loc k];
    end
end
P = length(pilot_loc);
Vj = sqrt(( 10 .^(SNR / 10) ) );
%% generation and test
QAMmod_order = 2;
serialpilot = randi([0,1],1, P * G * QAMmod_order);
parallelpilot = reshape(serialpilot, P * QAMmod_order, G);
parallelpilot_mod = qammod(parallelpilot,2^QAMmod_order,'InputType','bit','UnitAveragePower',true);
for g=1:G
    sequence_Phi_norm=parallelpilot_mod(:,g)/norm(parallelpilot_mod(:,g))*sqrt(P);
    sequence_Phi(:,g)=sequence_Phi_norm;
end
coefficients = 2 * randi([0 1], M, G) - 1;
coefficients=coefficients/sqrt(M);

n = [0:N-1]; 
kk = [0:N-1]; 
Wn = exp(-1i*2*pi/N); 
nk = n'*kk; 
Wnnk = Wn.^nk;
F = Wnnk ;
F = F(pilot_loc,1:L);
channel_H_delay = [0 randperm(L-1,4).*(1./B)].';
delay_H = round(channel_H_delay.*B) + 1;
channel_H = zeros(L,M);
theta=(sqrt(3)/2) * (2 * rand(length(delay_H),1) - 1);

for l = 1:length(delay_H)
    a_l = exp( -(0: M-1).'*1j*pi * theta(l) );
    h_l=exp( 1j*2*pi * rand() ) .* exp( -rand( ) );
    abs(h_l);
    g_1=h_l*a_l;
    channel_H(delay_H(l),:) =g_1.' ;
end

F_truncated = F;
for g=1:G
    R= diag(sequence_Phi(:,g)) * F * channel_H * coefficients(:,g);
    noise = sqrt(1/2) * (randn(P,1) + 1i * randn(P,1));
    y_received = R + 1 / Vj * noise;
    Y(:,g)=y_received;
end

%% received signals
H_true = channel_H ;




