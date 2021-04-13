function [check_flag, correct_flag, path_rank, dec_seq] = SLVD_dual_trellis(dual_trellis, dual_terminations, rxSig, crc_poly, tx_seq, Psi)

%
% This function is to perform soft serial list Viterbi decoding (SLVD) on the 
% dual trellis of a (n, n-1, v) convolutional encoder. The detailed dual trellis
% construction can be found in [1].
%
% References:
%   [1] T. Yamada et al., "A new maximum likelihood decoding of high rate
%       convolutional codes using a trellis", Elec. Commun. in Japan, 1983.
%
% Input parameters:
%   1) dual_trellis, a struct denoting the overall dual trellis
%   2) dual_terminations: a numState-by-numTermination matrix, with entry in
%       octal. Note that "numState" is the number of states in primal
%       trellis. The states in primal trellis are the first half states of 
%       dual trellis.
%   3) rxSig: a length n*h vector denoting the received real signals, where
%           h is the number of input-output streams.
%   4) crc_poly: CRC polynonmial in binary, degree from lowest to highest
%   5) tx_seq: a length (n-1)*h binary vector denoting the input sequence
%           to the (n, n-1, v) convolutional encoder. This is only to used
%           to determine "correct_flag".
%   6) Psi: a scalar denoting the constrained maximum list size.
%
% Output parameters:
%   1) check_flag: true if the decoded input sequence is divisible by
%           "crc_poly".
%   2) correct_flag: true if the decoded input sequence is indentical to
%           "tx_seq".
%   3) path_rank: a scalar denoting the terminating list rank at which
%           divisibility is first met.
%   4) dec_seq: the decoded input sequence that first passes CRC check.
%
% Remarks:
%   1) The decoder only works for (n, n-1, v) conv. encoder.
%   2) The function requires "PriorityQueue" class.
%
% Written by Hengjie Yang (hengjie.yang@ucla.edu)   04/08/21.
%



% Load important parameters
numStates = dual_trellis.numStates; % the # states on dual trellis
numTermination = size(dual_terminations, 2); % the # transitions for termination
n = dual_trellis.numOutputRails; % the # output rails of high-rate CC
rxSig = reshape(rxSig, n, []); % reshape the received vector
N = size(rxSig, 2); % the # transitions, each transition requires a dual trellis component

global Column
Column = cell(N+1, n);


% Re-organize termination symbols according to {iter, ii}{cur_state}
Termination_branches = cell(numTermination, n);

for init_state = 1:size(dual_terminations, 1)
    cur_state = init_state;
    for iter = 1:numTermination
        output_symbol = dual_terminations(init_state, iter);
        output_bits = dec2bin(oct2dec(output_symbol), n) - '0';
        for ii = 1:n
            if isempty(Termination_branches{iter, ii})
                Termination_branches{iter, ii} = cell(numStates, 1);
            end
            Termination_branches{iter, ii}{cur_state} = [Termination_branches{iter, ii}{cur_state}, output_bits(ii)]; 
            Termination_branches{iter, ii}{cur_state} = unique(Termination_branches{iter, ii}{cur_state}); % to make sure no duplicate bits are added
            next_state = dual_trellis.nextStates{ii}(cur_state, output_bits(ii)+1) + 1;
            cur_state = next_state;
        end
    end
end
            
             
            
            



% Step 1: Add-compare-select operations

% initialization
Column{1, 1} = cell(numStates, 1);
Column{1, 1}{1} = zeros(2, 2);
Column{1, 1}{1}(1,:) = [0, -1];
Column{1, 1}{1}(2,:) = [Inf, -1];

% In the following variables, the starting index is 1.

for iter = 1:N 
    for ii = 1:n
        Column{iter+floor(ii/n), mod(ii, n)+1} = cell(numStates, 1); % initialize the next depth
        
        if iter == 1 && ii == 1
            for input = 0:1
                next_state = dual_trellis.nextStates{ii}(1, input+1) + 1; 
                if next_state ~= 0
                    output_bit = input;
                    output_metric = -2*output_bit + 1; % BPSK rule
                    cum_metric = (rxSig(ii, iter)-output_metric)^2 + Column{iter, ii}{1}(1, 1); % cumulative metric
                    if isempty(Column{iter+floor((ii)/n), mod(ii, n)+1}{next_state})
                        Column{iter+floor(ii/n), mod(ii, n)+1}{next_state} = zeros(2, 2);
                        Column{iter+floor(ii/n), mod(ii, n)+1}{next_state}(1, :) = [cum_metric, 1]; % format: [metric, father_state]
                        Column{iter+floor(ii/n), mod(ii, n)+1}{next_state}(2, :) = [Inf, -1]; % format: [differential metric, father_state]
                    end
                end
            end
        elseif iter <= N - numTermination % can follow arbitrary branch
            for cur_state = 1:numStates
                if ~isempty(Column{iter, ii}{cur_state}) % a valid state
                    for input = 0:1
                        next_state = dual_trellis.nextStates{ii}(cur_state, input+1) + 1;
                        if next_state ~= 0
                            output_bit = input;
                            output_metric = -2*output_bit + 1; % BPSK rule
                            cum_metric = (rxSig(ii, iter)-output_metric)^2 + Column{iter, ii}{cur_state}(1, 1); % cumulative metric
                            if isempty(Column{iter+floor(ii/n), mod(ii, n)+1}{next_state})
                                Column{iter+floor(ii/n), mod(ii, n)+1}{next_state} = zeros(2, 2);
                                Column{iter+floor(ii/n), mod(ii, n)+1}{next_state}(1,:) = [cum_metric, cur_state];
                                Column{iter+floor(ii/n), mod(ii, n)+1}{next_state}(2,:) = [Inf, -1];
                            elseif Column{iter+floor(ii/n), mod(ii, n)+1}{next_state}(1,1) > cum_metric % previous path is worse
                                diff_metric = Column{iter+floor(ii/n), mod(ii, n)+1}{next_state}(1,1) - cum_metric;
                                state = Column{iter+floor(ii/n), mod(ii, n)+1}{next_state}(1,2);
                                Column{iter+floor(ii/n), mod(ii, n)+1}{next_state}(1,:) = [cum_metric, cur_state];
                                Column{iter+floor(ii/n), mod(ii, n)+1}{next_state}(2,:) = [diff_metric, state];
                            else % previous path remains optimal
                                diff_metric = cum_metric - Column{iter+floor(ii/n), mod(ii, n)+1}{next_state}(1,1);
                                Column{iter+floor(ii/n), mod(ii, n)+1}{next_state}(2,:) = [diff_metric, cur_state];
                            end
                        end
                    end
                end
            end
        else % enter the termination stage
            for cur_state = 1:numStates
                if ~isempty(Column{iter, ii}{cur_state}) && ~isempty(Termination_branches{iter-(N-numTermination), ii}{cur_state}) % a valid state
                    for jj = 1:size(Termination_branches{iter-(N-numTermination), ii}{cur_state}, 2)
                        input = Termination_branches{iter-(N-numTermination), ii}{cur_state}(jj); % only follow termination branches
                        next_state = dual_trellis.nextStates{ii}(cur_state, input+1) + 1;
                        if next_state ~= 0
                            output_bit = input;
                            output_metric = -2*output_bit + 1; % BPSK rule
                            cum_metric = (rxSig(ii, iter)-output_metric)^2 + Column{iter, ii}{cur_state}(1, 1); % cumulative metric
                            if isempty(Column{iter+floor(ii/n), mod(ii, n)+1}{next_state})
                                Column{iter+floor(ii/n), mod(ii, n)+1}{next_state} = zeros(2, 2);
                                Column{iter+floor(ii/n), mod(ii, n)+1}{next_state}(1,:) = [cum_metric, cur_state];
                                Column{iter+floor(ii/n), mod(ii, n)+1}{next_state}(2,:) = [Inf, -1];
                            elseif Column{iter+floor(ii/n), mod(ii, n)+1}{next_state}(1,1) > cum_metric % previous path is worse
                                diff_metric = Column{iter+floor(ii/n), mod(ii, n)+1}{next_state}(1,1) - cum_metric;
                                state = Column{iter+floor(ii/n), mod(ii, n)+1}{next_state}(1,2);
                                Column{iter+floor(ii/n), mod(ii, n)+1}{next_state}(1,:) = [cum_metric, cur_state];
                                Column{iter+floor(ii/n), mod(ii, n)+1}{next_state}(2,:) = [diff_metric, state];
                            else % previous path remains optimal
                                diff_metric = cum_metric - Column{iter+floor(ii/n), mod(ii, n)+1}{next_state}(1,1);
                                Column{iter+floor(ii/n), mod(ii, n)+1}{next_state}(2,:) = [diff_metric, cur_state];
                            end
                        end
                    end
                end
            end
        end
    end
end

% disp('Step 1 completed!');



% Step 2: the traceback of SLVD

decoded_input_sequences = [];
dec_seq = -1;
global Detour_tree
global Detour_array


Detour_tree = [];
node = struct('depth', N*n+5, 'father', -1);
Detour_tree = [Detour_tree; node];
Detour_array = [];

q = PriorityQueue([1, -1]);
counter = 0; % count the # tracebacks performed
path_rank = Psi;
check_flag = -1;
correct_flag = -1;

while counter < Psi
    counter = counter + 1;
    node = q.remove();
    id = node(1);
    [decoded_input_seq, next_GMD, next_state_ids] = trace_back_search(dual_trellis, N, id); % decoded_input_seq(1) first enters the conv. encoder
    decoded_input_sequences = [decoded_input_sequences; decoded_input_seq];
    if next_state_ids(1) ~= -1
        for ii = 1:size(next_state_ids, 2)
            q.insert([next_state_ids(ii), next_GMD]);
        end
    end
    
    % perform CRC verification
    temp = decoded_input_seq(1:(N-numTermination)*(n-1)); % extract the part without termination
    [~, remd] = gfdeconv(fliplr(temp), crc_poly);
    if any(remd) == 0
        check_flag = 1;
        path_rank = counter;
        dec_seq = decoded_input_seq;
        
        if all(decoded_input_seq == tx_seq)
            correct_flag = 1;
        else
            correct_flag = 0;
        end
        return
    end
end

check_flag = 0; % did not find any sequence that passes the CRC check within Psi trials
end





% Subfunction 1: trace back to identify the next most likely path on trellis
function [decoded_input_seq, next_GMD_value, next_path_ids] = trace_back_search(dual_trellis, N, id)

% The goal of this function is to identify the 'id'-th most likely path and
% the induced sequentially ordered paths from the 'id'-th path.
% 
% Input parameters:
%   1) dual_trellis: a struct as introduced before
%   2) N: the total # dual trellis components
%   3) id: the 'id'-th most likely path to be traced back
% Output parameters:
%   1) decoded_input_seq: the input sequence associated with the 'id'-th
%       most likely path.
%   2) GMD: the global merge difference (GMD) between the 'id'-th path and 1st
%       most likely path. This evalue is simply the metric difference
%       between these two paths.
%   3) next_path_ids: the sequence of indices representing the ranks of 
%       sequentially ordered most likely path induced from the 'id'-th path.
%

global Column
global Detour_array
decoded_input_seq = [];
n = dual_trellis.numOutputRails; % the length of the dual trellis component
decoded_codeword = -1*ones(1, N*n);
detour_trellis_depths = DFS(id); % to obtain the sequence of detour depths on the overall dual trellis
GMD_increment = 0; % a scalar denoting the GMD increment value on the go.
LMDs = zeros(1, N*n+1); % the LMD values along the path. For brevity, the size is equal to that of "Column".
LMDs(1) = Inf; % invalidate the 1st position as it is inaccessible.
cur_state = 1; % here, we automatically incremented to the true state.

for depth = N*n+1:-1:2
    iter = ceil(depth / n); % automatically incremented
    ii = depth - (iter - 1)*n; % automatically incremented
    LMDs(depth) = Column{iter, ii}{cur_state}(2, 1) + GMD_increment; % LMDs(i) represents the diverging point at i+1.
    if any(detour_trellis_depths == depth) % if depth is in the list, make a detour to the suboptimal path
        pre_state = Column{iter, ii}{cur_state}(2, 2); % automatically incremented
        GMD_increment = GMD_increment + Column{iter, ii}{cur_state}(2, 1);
    else % if not in the list, keep following the optimal path
        pre_state = Column{iter, ii}{cur_state}(1, 2);
    end

    % collect input bit on the go
    for input = 1:2 % automatically incremented
        prev_ii = mod(ii-1, n)+n*(1-sign(ii-1));
        if dual_trellis.nextStates{prev_ii}(pre_state, input) + 1 == cur_state % mod(ii-1, n)+n*(1-sign(ii-1)) is a fancy way to circularly decrement ii in the range [1, n]
            decoded_codeword(depth-1) = input - 1;
            cur_state = pre_state;
            break % stops as soon as the input bit leading to 'cur_state' is determined
        end
    end
end

% Step 2: Extract the information sequence from the codeword
decoded_codeword = reshape(decoded_codeword, n, []);
decoded_input_seq = decoded_codeword(2:end, :); % extract the systematic part
decoded_input_seq = reshape(decoded_input_seq, 1, []);



% Step 3: identify the next minimum GMD value to of the subsequent rank
% ordered paths. It is possible that multiple diverging spots achieve the
% same next minimum GMD value.

next_path_ids = [];

min_depth = min(detour_trellis_depths);
LMDs(min_depth:end) = Inf; % invalidate the impossible segment, since the values in this range has been used for previous paths.
Detour_array = [Detour_array; LMDs];

next_GMD_value = min(min(Detour_array)); % determine the next GMD value
% CAUTION: we also need to invalidate the value at which no divering point
% exists.

if next_GMD_value == Inf % did not find the valid next GMD value
    next_path_ids = -1;
    return
end

[rows, cols] = find(Detour_array == next_GMD_value); % rows represent the father path index at which the next rank ordered path is derived
Detour_array(rows, cols) = Inf; % replace with Inf to avoid finding the same value next time

next_path_ids = Add_nodes(rows, cols); % add newly found indices of next rank ordered oath to "Detour_tree"

end




function detour_path = DFS(path_id)

% This function is to read the sequence of detour depths that leads to
% "path_id"-th path.


global Detour_tree

detour_path = [];
cur = path_id;
if cur == 1
    return 
end

while cur ~= 1
    detour_path = [detour_path, Detour_tree(cur).depth];
    cur = Detour_tree(cur).father;
end

detour_path = detour_path(end:-1:1);


end




function [next_path_ids] = Add_nodes(father_path_ids, detour_depths)

% Add the newly found indices of next rank ordered path to "Detour_tree".
% Here, "father_path_ids" and "detour_depths" have the same vector length.

global Detour_tree

next_path_ids = [];
for ii = 1:size(father_path_ids, 1)
    node = struct('depth', detour_depths(ii), 'father', father_path_ids(ii));
    Detour_tree = [Detour_tree; node];
    next_path_ids = [next_path_ids, size(Detour_tree, 1)];
end


end
            
                        
                        
            
                
                
                
        





















