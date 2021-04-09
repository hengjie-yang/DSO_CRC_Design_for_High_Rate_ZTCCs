% This script is to verify if a randomly generated high-rate ZTCC can go
% from all-zero state to all-zero state on the overall dual trellis. If
% this is true for all generated codewords, then the dual trellis is
% correct.
%
% The script requires generated primal trellis and dual trellis.
%
% Test results:
%   04-08-21: all instances passed!
%
% Written by Hengjie Yang (hengjie.yang@ucla.edu)   04/08/21
%


clear all;
clc;


% Load primal trellis and dual trellis
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
load([path, fileName, '.mat'], 'myTrellis', 'Terminations');


fileName = ['Dual_trellis_v_',num2str(v), '_num_', num_string, 'den_', num2str(denominator)];
if ~exist([path, fileName, '.mat'], 'file')
    disp('File ',fileName, ' does not exist!');
    return
end
load([path, fileName, '.mat'], 'dual_trellis');



% System parameters
K = 6000; % information length

num_experiment = 1e3;
numTermination = size(Terminations, 2);

for iter = 1:num_experiment
    info_sequence = randi([0 1], 1, K);
    
    % add termination
    [~, last_state] = convenc(info_sequence, myTrellis);  
    for ii = 1:numTermination
        termination_bits = Terminations(last_state+1, ii);
        termination_bits = dec2bin(oct2dec(termination_bits), k) - '0';
        info_sequence = [info_sequence, termination_bits];
    end
    
    % convolutionally encode the CRC-coded sequence
    [codeword, fstate] = convenc(info_sequence, myTrellis);
    if fstate ~= 0
        error('ERROR: Incorrect termination!');
    end
      
    % walk along the overall dual trellis
    [fstate] = walk_dual_trellis(dual_trellis, codeword);
    
    if fstate ~= 0
        disp('Incorrect dual trellis!');
        return
    else
        disp(['Case: ',num2str(iter),' Instance passed.']);
    end
end




function [fstate] = walk_dual_trellis(dual_trellis, codeword)

cur_state = 0;
n = dual_trellis.numOutputRails;
codeword = reshape(codeword, n, []);
N = size(codeword, 2);

for iter = 1:N
    for ii = 0:n-1
        next_state = dual_trellis.nextStates{ii+1}(cur_state+1, codeword(ii+1, iter)+1);
        cur_state = next_state;
    end
end

fstate = cur_state;
end
        
    









