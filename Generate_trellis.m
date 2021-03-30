function [myTrellis] = Generate_trellis(v, numerators, denominator)

% The function generates the struct similar to 'poly2trellis' for 
% an (n, n-1, v-1) systematic conv. encoder. Namely,
% it includes the following fields:
%   1) numInputSymbols
%   2) numOutputSymbols
%   3) numStates
%   4) nextStates
%   5) outputs
%
% Input parameters of (n, n-1, v-1) systematic conv. encoder:
%   1) v-1: the overall constraint length
%   2) numerators: conventional octal form of numerators in n-th column
%   3) denominator: conventional octal form of den. in n-th column 
%
% Written by Hengjie Yang (hengjie.yang@ucla.edu)   03/28/21
%


k = length(numerators);
n = k + 1;

numInputSymbols = 2^k;
numOutputSymbols = 2^n;
numStates = 2^(v - 1);

nextStates = zeros(numStates, numInputSymbols);
outputs = zeros(numStates, numInputSymbols);

% Compute the 'nextStates' matrix

gs = zeros(k, v);

for iter=1:k
    gs(iter,:) = dec2bin(base2dec(num2str(numerators(iter)), 8)) - '0';
end
b = dec2bin(base2dec(num2str(denominator), 8)) - '0';


for currentState = 0:numStates-1
    cur_sigmas = dec2bin(currentState, v-1) - '0'; % from highest to lowest: [sigma_1, sigma_2, ..., sigma_v]
    for input = 0:numInputSymbols-1
        us = dec2bin(input, k) - '0'; % from highest to lowest: [u_1, u_2, ..., u_k]
        us = us'; % convert to a k-by-1 vector
        total_numerator = 0;
        for ii = 1:length(us)
            temp = gfconv(us(ii), fliplr(gs(ii, :))); % flip to lowest to highest
            total_numerator = gfadd(total_numerator, temp);
        end
        total_numerator = gfadd(total_numerator, [0, cur_sigmas]);
        [q, remd] = gfdeconv(total_numerator, fliplr(b)); % flip to lowest to highest
        next_sigmas = gfadd(zeros(1, v-1), remd);
        nextState = bi2de(fliplr(next_sigmas)); % due to the feature of 'bi2de', we need an extra fliplr. 
        % nextState represented in decimal from 0 to 2^v-1;
        nextStates(currentState+1, input+1) = nextState;
        lastbit = q;
        nextOutput = [us', lastbit];
        outputs(currentState+1, input+1) = str2double(dec2base(bi2de(fliplr(nextOutput)), 8));
    end
end

myTrellis.numInputSymbols = numInputSymbols;
myTrellis.numOutputSymbols = numOutputSymbols;
myTrellis.numStates = numStates;
myTrellis.nextStates = nextStates;
myTrellis.outputs = outputs;


% Search for shortest termination symbols
% here, we use backward breath first search to identify the termination
% symbols for each state.

num_transitions = ceil((v-1)/k);
Tree = zeros(numStates, 1); % Tree(i) = j means j is the father of i.
queue = [];
queue = [queue, 0];
vis = zeros(numStates, 1);
vis(1) = 1;

while ~isempty(queue)
    target_state = queue(1);
    queue = queue(2:end);
    for pre_state = 0:numStates-1
        if any(ismember(nextStates(pre_state+1, :), target_state))==1 && vis(pre_state+1) == 0
            vis(pre_state+1) = 1;
            queue = [queue, pre_state];
            Tree(pre_state+1) = target_state; 
        end
    end
end


Terminations = cell(numStates, 1);
Terminations{1} = [0];
for cur_state = 1:numStates-1
    cur = cur_state;
    while cur ~= 0
        fa_state = Tree(cur+1);
        index = find(nextStates(cur+1, :) == fa_state); 
        Terminations{cur_state+1} = [Terminations{cur_state+1}, ...
            str2double(dec2base(index-1, 8))]; % in octal
        cur = fa_state;
    end
end

% pad with zeros
for cur_state = 0:numStates-1
    len = length(Terminations{cur_state+1});
    if len < num_transitions
        Terminations{cur_state+1} = [Terminations{cur_state+1}, zeros(1, num_transitions-len)];
    end
end

myTrellis.numTermination = num_transitions;
myTrellis.terminations = Terminations;



% save results
path = './Simulation_results/';

num_string = '';
for iter = 1:k
    num_string = [num_string, num2str(numerators(iter)), '_'];
end

fileName = ['Trellis_v_',num2str(v), '_num_', num_string, 'den_', num2str(denominator)];
save([path, fileName, '.mat'], 'myTrellis');

            
    




