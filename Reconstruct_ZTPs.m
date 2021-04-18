function ZTP_node = Reconstruct_ZTPs(v, numerators, denominator, d_tilde, N)


% This function is to reconstruct all zero-terminated paths (ZTPs) from the
% irreducible error events (IEEs). Traditional method by Lou et al. does
% not adapt to the high-rate ZTCCs due to nontrivial terminations.
%
% Input parameters:
%   1) v: (v-1) denotes # memory elements in systematic feedback encoder.
%   2) numerators: a 1-by-k row vector with entries in octal, where k
%   denotes the # input rails
%   3) denominator: a scalar in octal
%   4) d_tilde: a scalar denoting the distance threshold (achievable)
%   5) N: a scalar denoting the primal trellis length
%
% Output parameters: ZTP_node a struct composed of following fields
%   1) list: a d_tilde-by-1 column vector denoting the list of length-kN
%       ZTPs arranged in ascending distances. The true distance is one less
%       than its index.
%   2) aggregate: a scalar denoting the number of length-kN ZTPs of
%       distance less than 'd_tilde'.
%
% Remarks:
%   1) Need to run "find_irreducible_error_event.m" if IEEs are not
%       generated before.
%   2) Need to run "Compute_ZTCC_weight_spectrum.m" if weight_node is not
%       generated before.
%   3) The distance index is true distance plus one
%
% Written by Hengjie Yang (hengjie.yang@ucla.edu)   04/17/21
%


% Step 1: load files
ZTP_node = {};

path = './Simulation_results/';
k = length(numerators);
num_string = '';
for iter = 1:k
    num_string = [num_string, num2str(numerators(iter)), '_'];
end

fileName = ['error_events_v_',num2str(v),'_num_',num_string,'den_',...
    num2str(denominator),'_d_tilde_',num2str(d_tilde)];

if ~exist([path, fileName, '.mat'], 'file')
    disp(['Error: the file ',fileName, ' does not exist!']);
    return
end

load([path, fileName, '.mat'], 'error_events', 'error_event_lengths');

fileName = ['weight_spectrum_v_',num2str(v), '_num_', num_string,'den_',...
    num2str(denominator),'_N_',num2str(N)];

if ~exist([path, fileName, '.mat'], 'file')
    disp(['Error: the file ',fileName, ' does not exist!']);
    return
end

load([path, fileName, '.mat'], 'weight_node');


full_weight_spectrum = weight_node.weight_spectrum;
d_max = size(full_weight_spectrum, 1);
if d_tilde > d_max - 1
    msg = 'Error: d_tilde exceeds the maximum possible distance!';
    disp(msg);
    return
end

 



% Step 2: use dynamic programming to reconstruct the length-kN ZTPs.
disp('Step 2: use dynamic programming to reconstruct ZTPs with objective length.');
ZTPs = cell(d_tilde+1, 1);

% preprocessing: need to add zeros(1,k) to IEEs
IEEs = cell(d_tilde+1, 1);
IEE_lens = cell(d_tilde+1, 1);
IEEs{1} = zeros(1, k);
IEE_lens{1} = k;
for dist = 1:d_tilde % true distance
    IEEs{dist+1} = error_events{dist};
    IEE_lens{dist+1} = error_event_lengths{dist};
end

Temp_ZTPs = cell(d_tilde+1, N+1);


for dist = 0:d_tilde % enumerate true objective distance
    disp(['     Current distance: ',num2str(dist)]);
    for len = 1:N % enumerate true objective trellis length
        for weight = dist:-1:0 % enumerate true IEE weight
            for ii = 1:size(IEEs{weight+1}, 1)
                l = IEE_lens{weight+1}(ii) / k; % convert to trellis length
                if weight == dist && l == len
                    Temp_ZTPs{dist+1, len+1} = int8(Temp_ZTPs{dist+1, len+1}); % save memory
                    Temp_ZTPs{dist+1, len+1} = [Temp_ZTPs{dist+1, len+1}; IEEs{weight+1}(ii, 1:(k*l))];
                elseif l < len && ~isempty(Temp_ZTPs{dist-weight+1, len-l+1})
                    [row, ~] = size(Temp_ZTPs{dist-weight+1, len-l+1});
                    Added_bits = repmat(IEEs{weight+1}(ii, 1:(k*l)), row, 1);
                    New_ZTPs = [Temp_ZTPs{dist-weight+1, len-l+1}, Added_bits];
                    Temp_ZTPs{dist+1, len+1} = int8(Temp_ZTPs{dist+1, len+1}); % save memory
                    Temp_ZTPs{dist+1, len+1} = [Temp_ZTPs{dist+1, len+1}; New_ZTPs];
                end
            end
        end
    end
end


% After building, extract the objective length ZTPs
for dist = 0:d_tilde
    if ~isempty(Temp_ZTPs{dist+1, N+1})
        ZTPs{dist+1} = Temp_ZTPs{dist+1, N+1};
    end
end

clearvars Temp_ZTPs


% Step 3: check if shifting is required.
need_shift_flag = 0;

for dist = 0:d_tilde
    if size(ZTPs{dist+1}, 1) ~= full_weight_spectrum(dist+1)
        need_shift_flag = 1;
        break
    end
end


aggregate = 0;
for dist = 0:d_tilde
    aggregate = aggregate + size(ZTPs{dist+1}, 1);
end

ZTP_node.list = ZTPs;
ZTP_node.aggregate = aggregate;


if need_shift_flag == 0
    disp('Congratulations! No need to shift before saving results!');
else
    disp('Sorry, the shifting operation is required...');
end

fileName = ['ZTP_node_v_',num2str(v), '_num_', num_string,'den_',...
    num2str(denominator),'_N_',num2str(N)];
save([path, fileName],'ZTP_node','-v7.3');



end