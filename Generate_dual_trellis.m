function [dual_trellis] = Generate_dual_trellis(v, numerators, denominator)

% This function is used to generate the dual trellis component of a given
% (n, n-1, v-1) conv. encoder, according to the method presented in [1]. The
% output is a struct with fields similar to that generated by
% "poly2trellis" function.
%
% References:
%   [1] T. Yamada et al., "A new maximum likelihood decoding of high rate
%       convolutional codes using a trellis", Elec. Commun. in Japan, 1983.
%   [2] S. Lin and D. J. Costello, "Error control coding", Pearson Education,
%       Inc. 2004.
%
% Input parameters:
%   1) v-1: the overall constraint length
%   2) numerators: notation consistent with [Lin and Costello, p. 539].
%       Namely, octals represent [h^(n-1)(D), ..., h^(1)(D)].
%   3) denominator: notation consistent with [Lin and Costello, p. 539].
%       Namely, octal represents h^(0)(D).
%
% Output parameters: dual_trellis, a struct composed of fields below
%   1) numOutputRails: the dual trellis length, i.e., the number of output rails
%   2) numInputSymbols: the maximum number of allowable inputs per state,
%       the input at depth j corresponds to output bit y_k^{(j)}.
%   3) numStates: the maximum number of allowable starting states per
%       depth, which is equal to 2^v. Note that in the (n-1)-th transition, 
%       s^(n-1) connects to s^(0).
%   4) nextStates: a numOutRails-by-1 cell, each cell is a numStates-by-numInputSymbols
%       matrix, with each entry denoting the index of the next state (starting from 0).
%       Prohibited start state will have its row as -2.
%       Inaccessible next state from a valid start state is marked as -1.
%
% Remarks:
%   1) Due to the construction method, the output bit of a branch is equal
%       to the input bit. Hence, there is no need to preserve the "outputs".
%   2) "nextStates" has the wrap-around feature at j=n-1 in that the
%       transition is from state s^(n-1) to s^(0).
%   3) The maximum instant response factor 'c' will be computed from
%       "numerators".
%   4) Format of the state index: Let s^(j) = (s_{v-1}^(j), s_{v-2}^(j),...,
%       s_0^(j)), j = 0,1,...,n-1, be the state vector representing the 
%       partial sums of (v+1) adders in the observer canonical form of H(D). 
%       Then its decimal representation plus one is the state index. Note that 
%       (v-1) is the overall constraint length.
%
% Written by Hengjie yang (hengjie.yang@ucla.edu)   04/08/21
%

dual_trellis = [];
numOutputRails = length(numerators) + 1;
numInputSymbols = 2;
numStates = 2^v; % Note that (v-1) represents the overall constraint length.

nextStates = cell(numOutputRails, 1);

for iter = 1:numOutputRails
    nextStates{iter} = zeros(numStates, numInputSymbols);
end


% Step 1: Re-arrange, expand octals into binary vector
hs = zeros(numOutputRails, v); % degree from highest to lowest. the i-th row corresponds to h^(i)(D)
for iter = 1:numOutputRails
    if iter == 1
        hs(iter, :) = dec2bin(base2dec(num2str(denominator), 8), v) - '0';
    else
        hs(iter, :) = dec2bin(base2dec(num2str(numerators(numOutputRails+1-iter)), 8), v) - '0';
    end
end
disp('Step 1 completed');


% Step 2: Compute the maximum instant response factor
index = find(hs(:, v)==1);
c = index(end) - 1; % typically a number greater than 0
disp('Step 2 completed');


num_half_states = numStates / 2;
cur_valid_states = zeros(numStates, 1);
next_valid_states = zeros(numStates, 1);

for ii = 0:num_half_states-1
    cur_valid_states(ii+1) = 1;
end

% Step 3: Compute the "nextStates" matrix
for j = 0:numOutputRails-1 
    if j ~= c && j ~= numOutputRails-1 % both branches available, no need to right shift
        for cur_state = 0:numStates - 1
            if cur_valid_states(cur_state+1) == 1
                cur_state_binary = dec2bin(cur_state, v) - '0';
                for input = 0:numInputSymbols - 1
                    next_state_binary = gfadd(cur_state_binary, input*hs(j+1, :));
                    next_state = bin2dec(num2str(next_state_binary));
                    next_valid_states(next_state + 1) = 1;
                    nextStates{j+1}(cur_state+1, input+1) = next_state;
                end
            else
                nextStates{j+1}(cur_state+1, :) = -2*ones(1, numInputSymbols);
            end
        end
    elseif j == c && j ~= numOutputRails-1 % single branch available, no need to right shift
        for cur_state = 0:numStates-1
            if cur_valid_states(cur_state+1) == 1
                cur_state_binary = dec2bin(cur_state, v) - '0';
                next_state_binary = gfadd(cur_state_binary, cur_state_binary(end)*hs(j+1, :));
                next_state = bin2dec(num2str(next_state_binary));
                next_valid_states(next_state + 1) = 1;
                nextStates{j+1}(cur_state+1, cur_state_binary(end)+1) = next_state;
                nextStates{j+1}(cur_state+1, 2-cur_state_binary(end)) = -1; % the other input is automatic inaccessible
            else
                nextStates{j+1}(cur_state+1, :) = -2*ones(1, numInputSymbols);
            end
        end
    elseif j ~= c && j == numOutputRails-1 % both branches available, need to right shift
        for cur_state = 0:numStates-1
            if cur_valid_states(cur_state+1) == 1
                cur_state_binary = dec2bin(cur_state, v) - '0';
                for input = 0:numInputSymbols - 1
                    next_state_binary = gfadd(cur_state_binary, input*hs(j+1, :));
                    next_state_binary = [0, next_state_binary(1:end-1)]; % perform right shifting
                    next_state = bin2dec(num2str(next_state_binary));
                    next_valid_states(next_state + 1) = 1;
                    nextStates{j+1}(cur_state+1, input+1) = next_state;
                end
            else
                nextStates{j+1}(cur_state+1, :) = -2*ones(1, numInputSymbols);
            end
        end
    else % j == c && j == numOutputRails-1, single branch available, need to right shift
        for cur_state = 0:numStates-1
            if cur_valid_states(cur_state+1) == 1
                cur_state_binary = dec2bin(cur_state, v) - '0';
                next_state_binary = gfadd(cur_state_binary, cur_state_binary(end)*hs(j+1, :));
                next_state_binary = [0, next_state_binary(1:end-1)]; % perform right shifting
                next_state = bin2dec(num2str(next_state_binary));
                next_valid_states(next_state + 1) = 1;
                nextStates{j+1}(cur_state+1, cur_state_binary(end)+1) = next_state;
                nextStates{j+1}(cur_state+1, 2-cur_state_binary(end)) = -1; % the other input is automatic inaccessible
            else
                nextStates{j+1}(cur_state+1, :) = -2*ones(1, numInputSymbols);
            end
        end
    end
    
    % reset "cur_valid_states" and "next_valid_states"
    cur_valid_states = next_valid_states;
    next_valid_states = zeros(numStates, 1);
end



% Step 4: save results
dual_trellis.numOutputRails = numOutputRails;
dual_trellis.numInputSymbols = numInputSymbols;
dual_trellis.numStates = numStates;
dual_trellis.nextStates = nextStates;

path = './Simulation_results/';
num_string = '';
for iter = 1:numOutputRails-1
    num_string = [num_string, num2str(numerators(iter)), '_'];
end

fileName = ['Dual_trellis_v_',num2str(v), '_num_', num_string, 'den_', num2str(denominator)];
save([path, fileName, '.mat'], 'dual_trellis');

end

        
        
        
            
        



    






