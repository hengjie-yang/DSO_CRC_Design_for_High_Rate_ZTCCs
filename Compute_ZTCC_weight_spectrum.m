function weight_node = Compute_ZTCC_weight_spectrum(v, numerators, denominator, N)

% 
%   This function computes the weight spectrum of a given high-rate ZTCC
%   of length N
%
%   Inputs:
%       1) v-1: the overall constraint length
%       2) numerators: conventional octal form of numerators in 1st column
%       3) denominator: conventional octal form of den. in 1st column
%       4) N: the trellis length
%
%   Outputs: weight_node, a struct that includes
%       1) weight_spectrum: a (d_max+1)-by-1 column vector denoting the #
%           codewords of weight i. Index 'i' represents weight 'i-1'.
%       2) overall_weight_function: a polynomial representation of WEF
% 
%   Written by Hengjie Yang (hengjie.yang@ucla.edu) 03/30/21
%


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


NumStates = myTrellis.numStates;
NumInputSymbols = myTrellis.numInputSymbols;
T = zeros(NumStates, NumStates);
T = sym(T);
syms X;


% Step 1: compute the one-step transfer function
disp('Step 1: Compute the one-step transfer function.');
for cur_state = 1:NumStates
    for input = 1:NumInputSymbols
        next_state = myTrellis.nextStates(cur_state,input) + 1;
        output_weight = sum(dec2bin(oct2dec(myTrellis.outputs(cur_state, input)))-'0');     
        output_symbol = convert_to_symbol(output_weight);
        T(cur_state, next_state) = output_symbol;
    end
end

disp(T);

% Step 2: Compute the weight enumerating function for finite-length TBCC
disp('Step 2: compute the weight enumerating function for each starting state.');
B = eye(NumStates);
B = sym(B);

for iter = 1:N
    disp(['Current depths: ',num2str(iter)]);
    B = B*T;
    B = expand(B);
end


% Step 3: Compute the overall weight enumerating function for finite-length
% TBCC
disp('Step 3: Compute the overall weight enumerating function.');

Poly = B(1, 1);
disp(Poly);
weight_spectrum = coeffs(Poly,'All');
weight_spectrum = fliplr(weight_spectrum);
weight_spectrum = double(weight_spectrum);
weight_spectrum = weight_spectrum';

weight_node.weight_spectrum = weight_spectrum;
weight_node.overall_weight_function = Poly;


% Save results

fileName = ['weight_spectrum_v_',num2str(v), '_num_', num_string,...
    'den_', num2str(denominator),'_N_',num2str(N),'.mat'];

save([path, fileName],'weight_node','-v7.3');




function chr = convert_to_symbol(weight)

chr = 0; %invalid string
syms X

if weight == 0
    chr = 1;
elseif weight == 1
    chr = X;
else
    chr = X^weight;
end
    

