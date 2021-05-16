function [Poly_node] = Search_DSO_CRC_poly(v, numerators, denominator, d_tilde, N, m, base)


% This function is to identity the DSO CRC polynomial that maximizes the
% minimum distance for the given high-rate ZTCC.
%
% Input parameters:
%   1) v: (v-1) denotes the # memory elements
%   2) numerators: a 1-by-k row vector with entries in octal, where k
%   denotes the # input rails
%   3) denominator: a scalar in octal
%   4) d_tilde: a scalar denoting the true distance threshold
%   5) N: a scalar denoting the primal trellis length
%   6) m: the objective CRC degree
%   7) base: the representation of CRC polynomial; default value is 16
%
% Output parameters: Poly_node, a struct including the following fields
%   1) success_flag: true if the unique DSO CRC polynomial is identified.
%   2) crc_gen_polys: the final list of CRC candidates; the 
%                list size = 1 (i.e., DSO CRC poly) if success_flag = True, 
%                and list size >1 otherwise;
%   3) stopped_distance: -1 if success_flag = True, and a scalar
%                indicating the stopped distance at which the DSO CRC is
%                still not finalized if success_flag = False;
%   4) crc_distance: a positive integer if undetected errors are
%                found within "d_tilde", and -1 otherise;
%
% Written by Hengjie Yang (hengjie.yang@ucla.edu)   04/18/21
%

tic

Poly_node = {};
stopped_distance = -1;
success_flag = false;
min_dist = -1;
k = length(numerators); % # input rails
mu = ceil((v-1)/k); % the # termination transitions

if nargin < 7
    base = 16;
end

% Step 1: load files
path = './Simulation_results/';
num_string = '';
for iter = 1:k
    num_string = [num_string, num2str(numerators(iter)), '_'];
end
fileName = ['ZTP_node_v_',num2str(v), '_num_', num_string,'den_',...
    num2str(denominator), '_d_',num2str(d_tilde),'_N_',num2str(N)];

if ~exist([path, fileName, '.mat'], 'file')
    disp(['Error: the file ',fileName, ' does not exist!']);
    return
end

load([path, fileName, '.mat'], 'ZTP_node');

ZTPs = ZTP_node.list; % the index is true distance plus one



% Step 2: search the DSO CRC polynomial

List_size = 2^(m-1);
Candidate_CRCs = dec2bin(0:List_size-1) - '0';
Candidate_CRCs = [ones(List_size,1), Candidate_CRCs, ones(2^(m-1),1)]; % degree order from highest to lowest
Undetected_spectrum = inf(List_size, 1); % each column represents the undected spectrum

Candidate_poly_octal=dec2base(bin2dec(num2str(Candidate_CRCs)),base); % octal form


if m == 1 % special case
    Candidate_CRCs = [1 1];
    Undetected_spectrum = inf(List_size, 1); % each column represents the undected spectrum
    Candidate_poly_octal=dec2base(bin2dec(num2str(Candidate_CRCs)),base); % octal form
end

mask = true(size(Candidate_CRCs,1),1);
locations = find(mask == true);

crc_gen_polys = [];
crc_gen_poly_vecs = [];


for dist = 2:d_tilde + 1 % incremented distance
    Undetected_spectrum = [Undetected_spectrum, inf(List_size, 1)];
    weight_vec = zeros(size(locations, 1), 1);
    if ~isempty(ZTPs{dist})
        parfor ii = 1:size(locations, 1)
            weight_vec(ii) = Check_divisible_by_distance(Candidate_CRCs(locations(ii),:), ZTPs{dist}, k, mu);
        end
        
        % copy the result to undetected_spectrum
        for ii = 1:size(locations, 1)
            Undetected_spectrum(locations(ii), dist) = weight_vec(ii);
        end
        
        min_weight = min(weight_vec);
        locations = locations(weight_vec == min_weight);
        disp(['    Current distance: ',num2str(dist-1),' number of candidates: ',...
            num2str(size(locations,1))]);
        if length(locations) == 1
            crc_gen_polys = Candidate_poly_octal(locations(1),:);
            crc_gen_poly_vecs = Candidate_CRCs(locations(1),:);
            success_flag = true;
            break
        end
    end
    
    if dist == d_tilde + 1 && length(locations) > 1
        crc_gen_polys = Candidate_poly_octal(locations,:);
        crc_gen_poly_vecs = Candidate_CRCs(locations,:);
        stopped_distance = d_tilde + 1; % the incremented distance
        disp(['    d_tilde is insufficient to find the DSO CRC... ']);
        disp(['    Stopped distance: ',num2str(stopped_distance),...
            ' # of candidate polynomials: ',num2str(size(crc_gen_polys,1))]);
    end
end


% Step 3: identify the minimum distance
if success_flag == true
    disp('Step 4: Identify the minimum undetected distance by the DSO CRC.');
    for dist = 2:d_tilde + 1 % the incremented distance 
        if ~isempty(ZTPs{dist})
            w = Check_divisible_by_distance(crc_gen_poly_vecs(1,:), ZTPs{dist}, k, mu);
            if w > 0
                min_dist = dist - 1;
                disp(['    DSO CRC polynomial: ',num2str(crc_gen_polys(1,:))]);
                disp(['    Minimum undetected distance: ',num2str(min_dist)]);
                break
            end
            if w == 0 && dist == d_tilde + 1
                disp('    d_tilde is insufficient to determine the minimum undetected distance.');
            end
        end
    end
end

% Save results
Poly_node.success_flag = success_flag;
Poly_node.crc_gen_polys = crc_gen_polys;
Poly_node.stopped_distance = stopped_distance;
Poly_node.crc_distance = min_dist;

file_name = ['Poly_node_high_rate_ZTCC_v_',num2str(v), '_num_', num_string,'den_',...
    num2str(denominator), '_d_',num2str(d_tilde),'_N_',num2str(N),'_m_',num2str(m),'.mat'];
save([path, file_name],'Poly_node','-v7.3');

timing = toc;
disp(['Execution time: ',num2str(timing),'s']);


end




function weight = Check_divisible_by_distance(poly_vec, error_events, k, mu)

% This function computes the undetected weight for "poly_vec" based on "error_events".
% Input parameters:
%   1) poly_vec: a binary vector with descending degree
%   2) error_events: the set of length-kN input sequences with the same output
%       distance
%   3) k: the # input rails
%   4) mu: the # termination transitions

weight = 0;
poly_vec = fliplr(poly_vec); % flip degree order from lowest to highest

for i = 1:size(error_events, 1)
    temp = double(error_events(i,:));
    temp = temp(1:end-k*mu); % remove the termination part
    temp = fliplr(temp); % flip input sequence degree to "lowest to highest"
    [~, remd] = gfdeconv(temp,poly_vec, 2);
    if remd == 0
        weight = weight + 1;
    end
end

end






