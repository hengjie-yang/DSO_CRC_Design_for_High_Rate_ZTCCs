function [error_events, error_event_lengths] = find_irreducible_error_event(v, numerators, denominator, d_tilde)

%
% The function aims at gathering all error events of the high-rate trellis
% up to distance 'd_tilde'.
%
% Inputs to specify an (n, n-1, v) systematic conv. encoder:
%   1) v: the overall constraint length + 1
%   2) numerators: conventional octal form of numerators in n-th column
%   3) denominator: conventional octal form of den. in n-th column 
%   
% Written by Hengjie Yang (hengjie.yang@ucla.edu)   03/28/21
%



% Step 1: load files

path = './Simulation_results/';
k = length(numerators);
num_string = '';
for iter = 1:k
    num_string = [num_string, num2str(numerators(iter)), '_'];
end

fileName = ['Trellis_v_',num2str(v), '_num_', num_string, 'den_', num2str(denominator)];

if ~exist([path, fileName, '.mat'], 'file')
    disp(['Error: the file ',fileName, ' does not exist!']);
    return
end

load([path, fileName, '.mat'], 'myTrellis');

numStates = myTrellis.numStates;
numInputSymbols = myTrellis.numInputSymbols;
numOutputSymbols = myTrellis.numOutputSymbols;
Terminations = myTrellis.terminations;



% Step 2: Viterbi search of irreducible error events

MaxIteration = 100;
Zero_state = cell(MaxIteration, 1);
Column = cell(2, 1);

for iter = 0:MaxIteration-1
    disp(['Current trellis depth: ', num2str(iter)]);
    Column{mod(iter+1, 2)+1} = cell(numStates, 1);
    if iter == 0
        for input = 1:numInputSymbols-1 % deviate from all-zero branch
            next_state = myTrellis.nextStates(1, input+1);
            if isempty(Column{mod(iter+1, 2)+1}{next_state+1})
                Column{mod(iter+1, 2)+1}{next_state+1} = cell(d_tilde, 1);
            end
            weight = sum(dec2bin(oct2dec(myTrellis.outputs(1,input+1)))-'0');
            Column{mod(iter+1, 2)+1}{next_state+1}{weight} = ...
                [Column{mod(iter+1, 2)+1}{next_state+1}{weight}; input];
        end
    else
        for cur_state = 1:numStates-1
            if ~isempty(Column{mod(iter, 2)+1}{cur_state+1})
                for dist = 1:d_tilde % since the path is nonzero, the inputs are also nonzero, true distance >= 1
                    if ~isempty(Column{mod(iter, 2)+1}{cur_state+1}{dist})
                        for input = 0:numInputSymbols-1 % free to go anywhere
                            next_state = myTrellis.nextStates(cur_state+1, input+1);
                            if isempty(Column{mod(iter+1, 2)+1}{next_state+1})
                                Column{mod(iter+1, 2)+1}{next_state+1} = cell(d_tilde, 1);
                            end
                            weight = sum(dec2bin(oct2dec(myTrellis.outputs(cur_state+1, input+1)))-'0');
                            temp = Column{mod(iter, 2)+1}{cur_state+1}{dist};
                            temp = [temp, input*ones(size(temp, 1), 1)];
                            if dist + weight <= d_tilde
                                Column{mod(iter+1, 2)+1}{next_state+1}{dist+weight} = ...
                                    [Column{mod(iter+1, 2)+1}{next_state+1}{dist+weight}; temp];
                            end
                        end
                    end
                end
            end
        end
    end
    Zero_state{iter+1} = Column{mod(iter+1, 2)+1}{1};
end

error_events = cell(d_tilde, 1);
error_event_lengths = cell(d_tilde, 1);

for iter = 1:MaxIteration
    if ~isempty(Zero_state{iter})
        for dist = 1:d_tilde
            if ~isempty(Zero_state{iter}{dist})
                if isempty(error_events{dist})
                    error_events{dist} = Zero_state{iter}{dist};
                else
                    error_events{dist}=[error_events{dist},...
                        zeros(size(error_events{dist},1), iter-size(error_events{dist},2))];
                    error_events{dist}=[error_events{dist};Zero_state{iter}{dist}];
                end
                error_event_lengths{dist}=[error_event_lengths{dist};...
                        iter*ones(size(Zero_state{iter}{dist},1),1)];
            end
        end
    end
end


% Step 3: Save results
fileName = ['error_events_v_',num2str(v),'_num_',num_string,'den_',num2str(denominator)];
save([path, fileName, '.mat'], 'error_events', 'error_event_lengths', '-v7.3');
                            
                        
            












