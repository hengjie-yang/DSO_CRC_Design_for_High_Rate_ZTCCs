% This script is to simulate the performance of the high-rate ZTCC decoded
% with brute-force serial list Viterbi decoding (i.e., ML decoding). The
% goal is to have a baseline result that can be used to compare with list
% decoding over the dual trellis.
%
% The output DistTable has the following format: for 1<=i<= Psi
%      1) (i, 1) denotes # total instances at list rank 'i'.
%      2) (i, 2) denotes # correct decoding at list rank 'i'.
%      3) (i, 3) denotes # undetected errors at list rank 'i'.
%   (Psi+1, 1) denotes # NACKs for list ranks [1, Psi].
%
% Written by Hengjie Yang (hengjie.yang@ucla.edu)   04/06/21
%

clear all;
clc;

% Load trellis
v = 7;
numerators = [107, 135, 133];
denominator = 141;
k = length(numerators);

num_string = '';
for iter = 1:k
    num_string = [num_string, num2str(numerators(iter)), '_'];
end

path = './Simulation_results/';
fileName = ['Trellis_v_',num2str(v), '_num_', num_string, 'den_', num2str(denominator)];

if ~exist([path, fileName, '.mat'], 'file')
    disp('File ',fileName, ' does not exist!');
    return
end
load([path, fileName, '.mat'], 'myTrellis', 'Terminations', 'Dual_terminations');

fileName = ['Dual_trellis_v_',num2str(v), '_num_', num_string, 'den_', num2str(denominator)];

if ~exist([path, fileName, '.mat'], 'file')
    disp('File ',fileName, ' does not exist!');
    return
end
load([path, fileName, '.mat'], 'dual_trellis');


% System parameters
K = 6; % information length
m = 6; % CRC degree
snr_dBs = 1;
crc_gen_poly = '6F'; % degree from highest to lowest
poly = dec2base(base2dec(crc_gen_poly, 16), 2) - '0';
poly = fliplr(poly); % degree from lowest to highest

P_UE = zeros(1, length(snr_dBs));


% Step 1: Generate all high-rate codewords for brute-force list decoding
K_high = K + m; % the information length of the high-rate code
num_codewords = 2^(K_high);
Psi = min(num_codewords, 1e5); % maximum list size
temp = 0:num_codewords-1;
numTermination = size(Terminations, 2);
N = (K_high + numTermination*k)/k*(k+1); % the blocklength

info_sequences = dec2bin(temp) - '0';
high_rate_codewords = zeros(num_codewords, N);

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
end
disp('Step 1 completed');



% Step 2: Simulation part
DistTable = cell(length(snr_dBs), 1);
Ave_list_ranks = zeros(length(snr_dBs), 1);
P_UEs = zeros(length(snr_dBs), 1);
P_NACKs = zeros(length(snr_dBs), 1);

for iter = 1:size(DistTable, 1)
    DistTable{iter} = zeros(Psi+1, 3);
end


for iter = 1:size(snr_dBs, 2)
    snr = 10^(snr_dBs(iter)/10);
    
    num_UE = 0;
    num_NACK = 0;
    num_trial = 0;
    while num_UE < 10 || num_trial < 1e4
        num_trial = num_trial + 1;
        info_sequence = randi([0 1], 1, K);
%         info_sequence = ones(1, K);
        
        % add CRC behind the message
        msg_temp = fliplr(info_sequence); % degree from low to high
        [~, remd] = gfdeconv([zeros(1, m), msg_temp], poly, 2); 
        crc_coded_sequence = gfadd(remd, [zeros(1, m), msg_temp]);
        crc_coded_sequence = fliplr(crc_coded_sequence); %degree from high to low
        
        % add termination
        [~, last_state] = convenc(crc_coded_sequence, myTrellis);  
        for ii = 1:numTermination
            termination_bits = Terminations(last_state+1, ii);
            termination_bits = dec2bin(oct2dec(termination_bits), k) - '0';
            crc_coded_sequence = [crc_coded_sequence, termination_bits];
        end
        
        % convolutionally encode the CRC-coded sequence
        [codeword, fstate] = convenc(crc_coded_sequence, myTrellis);
        if fstate ~= 0
            error('ERROR: Incorrect termination!');
        end
        
        % BPSK modulation
        txSig = -2*codeword + 1;
        
        % Send txSig over the AWGN channel
        rxSig = awgn(txSig, snr_dBs(iter), 'measured');
        
        % brute-force soft SLVD
%         [check_flag, correct_flag, path_rank] = ...
%                 SLVD_brute_force(rxSig, high_rate_codewords, Psi, crc_coded_sequence, poly);
%         disp(['Brute-force: SNR (dB): ', num2str(snr_dBs(iter)), ' # trials: ',num2str(num_trial),...
%             ' # errors: ', num2str(num_UE), ' check: ',num2str(check_flag)...
%             ' correct: ',num2str(correct_flag), ' list_rank: ', num2str(path_rank)]);
        
        [check_flag, correct_flag, path_rank, ~] = ...
            SLVD_dual_trellis(dual_trellis, Dual_terminations, rxSig, poly, crc_coded_sequence, Psi);
            
        disp(['Dual trellis: SNR (dB): ', num2str(snr_dBs(iter)), ' # trials: ',num2str(num_trial),...
            ' # errors: ', num2str(num_UE), ' check: ',num2str(check_flag)...
            ' correct: ',num2str(correct_flag), ' list_rank: ', num2str(path_rank)]);
        
        % keep record of the result
        if check_flag == 0
            num_NACK = num_NACK + 1;
        elseif check_flag == 1 && correct_flag == 0
            num_UE = num_UE + 1;
        end
        
        if check_flag == 1 && path_rank <= Psi % the 2nd condition is redundant but necessary
            if correct_flag == 1
                DistTable{iter}(path_rank, 2) = DistTable{iter}(path_rank, 2) + 1;
            else
                DistTable{iter}(path_rank, 3) = DistTable{iter}(path_rank, 3) + 1;
            end
            DistTable{iter}(path_rank, 1) = DistTable{iter}(path_rank, 1) + 1;
        else
            DistTable{iter}(Psi+1, 1) = DistTable{iter}(Psi+1, 1) + 1;
        end  
    end
end



% Step 3: Process data to compute the average P_UE, P_NACK and list rank
list_ranks = 1:Psi;

for ii = 1:size(snr_dBs, 2)
    tot = sum(DistTable{ii}(1:Psi+1, 1));
    overall_distribution = DistTable{ii}(1:Psi, 1);
    overall_distribution(Psi) = overall_distribution(Psi) + DistTable{ii}(Psi+1, 1);
    overall_distribution = overall_distribution / tot;
    overall_distribution = overall_distribution';
    Ave_list_ranks(ii) = sum(list_ranks.*overall_distribution);
    
    P_UEs(ii) = sum(DistTable{ii}(1:Psi, 3)) / tot;
    P_NACKs(ii) = DistTable{ii}(Psi+1, 1) / tot;
end


% Step 4: Save results
timestamp = datestr(now, 'mmddyy_HHMMSS');
path = './Simulation_results/';
save([path, timestamp, '_sim_data_CRC_ZTCC_v_',...
    num2str(v),'_num_', num_string, 'den_', num2str(denominator),'_CRC_',crc_gen_poly,...
    '_K_',num2str(K),'.mat'],'snr_dBs','DistTable','Ave_list_ranks','P_UEs','P_NACKs');








        
        
    









