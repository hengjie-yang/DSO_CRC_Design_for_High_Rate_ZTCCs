function [myTrellis, Terminations, Dual_terminations] = Generate_trellis(v, numerators, denominator)

% The function generates the struct similar to 'poly2trellis' for 
% an (n, n-1, v-1) systematic conv. encoder. Namely,
% the generated (primal) trellis includes the following fields:
%   1) numInputSymbols: 2^(n-1), where n = length(numerators) + 1;
%   2) numOutputSymbols: 2^n
%   3) numStates: 2^(v-1), where v is overall constraint length + 1;
%   4) nextStates: a numStates-by-numInputSymbols matrix;
%   5) outputs: a numStates-by-numInputSymbols matrix; each entry is in
%           octal.
%   6) Terminations: a numStates-by-numTermination matrix; each entry is in
%           octal.
%
% References:
%   [1] S. Lin and D. J. Costello, "Error control coding", Pearson Education,
%       Inc. 2004. 
%
% Format of the numerators: For (n, n-1, v) encoder, the notation is
% consistent with [Lin and Costello, page 540]. Namely,
%
%   numerators = [h^(n-1)(D), ..., h^(1)(D)].
%   denominator = h^(0)(D)
%
% Format of the state index: Let [s_{v-2}, s_{v-1}, ..., s_0] represents
%   the state vector in the observer canonical form of H(D). Then its
%   decimal representation plus one is the state index in MATLAB.
%
% This implies that c^(0) is the coded symbol and c^(1)(D),...,c^(n-1)(D)
% correspond to information sequences u^(0)(D), ..., u^(n-1)(D).
%
% The convolutional codeword is formed according to
%   [c^(0)_i, c^(1)_i, ..., c^(n-1)_i] at time instant i.
% 
%
% Input parameters of (n, n-1, v-1) systematic conv. encoder:
%   1) v-1: the overall constraint length
%   2) numerators: octal form (Lin and Costello) of numerators in n-th column
%   3) denominator: octal form (Lin and Costello) of den. in n-th column
%
% Output parameters:
%   1) myTrellis: a struct with fields same as output of 'poly2trellis'.
%   2) Terminations: a numState-by-numTermination matrix indicating the
%       sequential input termination symbols in octal.
%   3) Dual_terminations: a numState-by-numTermination matrix indicating the
%       sequential output termination symbols in octal.
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

gs = zeros(k, v); % degree from highest to lowest. the i-th row corresponds to u^(i)(D)

for iter=1:k
    gs(iter,:) = dec2bin(base2dec(num2str(numerators(n-iter)), 8), v) - '0';
end
b = dec2bin(base2dec(num2str(denominator), 8), v) - '0';


for currentState = 0:numStates-1
    cur_sigmas = dec2bin(currentState, v-1) - '0'; % from highest to lowest: [sigma_v, sigma_{v-1}, ..., sigma_1]
    for input = 0:numInputSymbols-1
        us = dec2bin(input, k) - '0'; % from highest to lowest: [u_1, u_2, ..., u_k]
        us = us'; % convert to a k-by-1 vector
        total_numerator = 0;
        for ii = 1:length(us)
            temp = gfconv(us(ii), gs(ii, :)); % flip to lowest to highest
            total_numerator = gfadd(total_numerator, temp);
        end
        total_numerator = gfadd(total_numerator, [0, cur_sigmas]);
        [q, remd] = gfdeconv(total_numerator, b); % flip to lowest to highest
        next_sigmas = gfadd(zeros(1, v-1), remd);
        nextState = bi2de(fliplr(next_sigmas)); % due to the feature of 'bi2de', we need an extra fliplr. 
        % nextState represented in decimal from 0 to 2^v-1;
        nextStates(currentState+1, input+1) = nextState;
        nextOutput = [q, us']; % the 1st output symbol is the coded bit
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

% change the format from a family of cells to a single matrix
temp = Terminations;
Terminations = zeros(numStates, num_transitions);
for iter = 1:numStates
    Terminations(iter, :) = temp{iter};
end


% Compute terminations for dual trellis
Dual_terminations = zeros(numStates, num_transitions);

for ii = 0:numStates-1
    cur_state = ii;
    for iter = 1:num_transitions
        input_symbol = Terminations(ii+1, iter); % in octal form
        input_bits = dec2bin(oct2dec(input_symbol), k) - '0';
        [output_bits, next_state] = convenc(input_bits, myTrellis, cur_state);
        output_symbol = str2double(dec2base(bi2de(fliplr(output_bits)), 8));
        Dual_terminations(ii+1, iter) = output_symbol;
        cur_state = next_state;
    end
end



% save results
path = './Simulation_results/';

num_string = '';
for iter = 1:k
    num_string = [num_string, num2str(numerators(iter)), '_'];
end

fileName = ['Trellis_v_',num2str(v), '_num_', num_string, 'den_', num2str(denominator)];
save([path, fileName, '.mat'], 'myTrellis', 'Terminations','Dual_terminations');

            
    




