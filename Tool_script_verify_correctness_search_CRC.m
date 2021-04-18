% Tool script: to verify the correctness of the "Search_DSO_CRC_poly.m"
% script. 
%
% Written by Hengjie Yang (hengjie.yang@ucla.edu)   04/06/21
%


% Here, we try to use brute force method to identify the optimal deg-6 CRC
% polynomial for K = 6, (4,3,6) ZTCC with H = [107, 135, 133, 141];

clear all;
clc;


set(0,'DefaultTextFontName','Times','DefaultTextFontSize',18,...
    'DefaultAxesFontName','Times','DefaultAxesFontSize',16,...
    'DefaultLineLineWidth',1,'DefaultLineMarkerSize',7.75);
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');



% Step 1: Generate all high-rate codewords (without CRC)

path = './Simulation_results/';
load([path, 'Trellis_v_7_num_107_135_133_den_141.mat'], 'myTrellis', 'Terminations');
load([path, 'RCU_and_MC_bound_N_28_K_9.mat'], 'gamma_s','rcu_bounds');
load([path, '040621_173925_sim_data_CRC_ZTCC_v_7_num_107_135_133_den_141_CRC_6F_K_6.mat'], 'snr_dBs','P_UEs');

k = 3; % the number of input rails
K = 6; % the information length of the low-rate code
m = 6; % the CRC degree
K_high = K + m; % the information length of the high-rate code

num_codewords = 2^(K_high);
temp = 0:num_codewords-1;
info_sequences = dec2bin(temp) - '0';
numTermination = size(Terminations, 2);
N = (K_high + numTermination*k)/k*(k+1); % the blocklength

high_rate_codewords = zeros(num_codewords, N);
High_rate_spectrum = zeros(N+1, 1);

for iter = 1:size(info_sequences, 1)
    info_sequence = info_sequences(iter, :);
    
    % add termination
    [~, last_state] = convenc(info_sequence, myTrellis);  
    for ii = 1:numTermination
        termination_bits = Terminations(last_state+1, ii);
        termination_bits = dec2bin(oct2dec(termination_bits), k) - '0';
        info_sequence = [info_sequence, termination_bits];
    end
    
    [codeword, fstate] = convenc(info_sequence, myTrellis);
    if fstate ~= 0
        error('ERROR: Incorrect termination!');
    end
    high_rate_codewords(iter, :) = codeword;
    weight = sum(codeword, 2);
    High_rate_spectrum(weight+1) = High_rate_spectrum(weight+1) + 1;
end
disp('Step 1 completed');



% Step 2: Use brute-force method to identify the optimal CRC polynomial
base = 16;
List_size = 2^(m-1);
Candidate_CRCs = dec2bin(0:List_size-1) - '0';
Candidate_CRCs = [ones(List_size,1), Candidate_CRCs, ones(2^(m-1),1)]; % degree order from highest to lowest
Candidate_poly_base = dec2base(bin2dec(num2str(Candidate_CRCs)),base);

Low_rate_spectra = zeros(List_size, N+1);
dmins = zeros(List_size, 1);
opt_dmin = -1;

for iter = 1:List_size
    Low_rate_spectra(iter, :)  = Compute_low_rate_spectra(Candidate_CRCs(iter, :),info_sequences, high_rate_codewords);
    d_min = find(Low_rate_spectra(iter, :) > 0);
    dmins(iter) = d_min(2) - 1;
    opt_dmin = max(opt_dmin, dmins(iter));
end
disp('Step 2 completed');

%% 
opt_index = -1;
index = find(dmins == opt_dmin);
if length(index) == 1
    disp(['Optimal CRC poly in hex: ', Candidate_poly_base(index, :),...
        ' minimum distance: ', num2str(opt_dmin)]);
    opt_index = index;
else
    opt_NNs = min(Low_rate_spectra(index, opt_dmin+1));
    index2 = find(Low_rate_spectra(index, opt_dmin+1) == opt_NNs);
    disp(['Optimal CRC poly in hex: ', Candidate_poly_base(index(index2), :),...
        ' minimum distance: ', num2str(opt_dmin)]);
    opt_index = index(index2);
end


% Step 3: Plot union bound of the optimal CRC-ZTCC code
SNRs = 1:0.5:4;
Union_bound = zeros(1, size(SNRs, 2));
if opt_index > 0
    for iter = 1:size(SNRs, 2)
        A = sqrt(10^(SNRs(iter)/10));
        dists = opt_dmin:N;
        spec = Low_rate_spectra(opt_index, opt_dmin+1:end);
        temp = spec.*qfunc(A.*sqrt(dists));
        Union_bound(iter) = sum(temp, 2);
    end
    
    figure;
    semilogy(SNRs, Union_bound, '-.'); hold on
    semilogy(snr_dBs, P_UEs, '-+'); hold on
    semilogy(gamma_s, rcu_bounds, '-'); hold on
    legend('Union bound', 'Simulated', 'RCU bound');
    grid on
    xlabel('SNR (dB)');
    ylabel('Probability of UE');
    title('K = 6, m = 6, (4, 3, 6) ZTCC');
end






function spec = Compute_low_rate_spectra(crc_poly, info_sequences, codewords)

% Input parameters:
%   1) crc_poly: from highest to lowest degree
%   2) info_sequence: index 1 corresponds to the first bit entering the conv.
%       encoder
%   3) k: the # input rails
%   4) numTermination: the # termination transitions
%
%

poly = fliplr(crc_poly); % degree from lowest to highest
N = length(codewords(1, :));
spec = zeros(1, N+1);

for iter = 1:size(info_sequences, 1)
    info_seq = info_sequences(iter, :);
    temp = fliplr(info_seq); % degree from lowest to highest
    [~, remd] = gfdeconv(temp, poly);
    if any(remd) == 0
        weight = sum(codewords(iter, :), 2);
        spec(weight+1) = spec(weight+1) + 1;
    end
end



end








