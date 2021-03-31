function Poly_node = Search_DSO_CRC_poly(v, numerators, denominator, d_tilde, K, m, base)

% This function is to identify the DSO CRC polynomial for the (n, n-1, v)
% systematic feedback conv. encoder that maximizes the minimum distance of 
% the resulting concatenated code. 
%
% Input parameters:
%   1) v-1: the overall constraint length
%   2) numerators: conventional octal form of numerators in n-th column
%   3) denominator: conventional octal form of den. in n-th column
%   4) K: the total information length of all rails.
%   5) m: the objective CRC polynomial degree
%   6) base: a scalar denoting the base of CRC presentation
%
% Remarks:
%   1) (K+m) must be divisible by (n-1). Otherwise, the function will
%   immediately return.
%   2) If the 1st condition is met, the trellis length N  = (k+m)/(n-1).
%   3) CRC is represented from highest to lowest degree.
%   4) The function requires output files from
%           find_irreducible_error_event.m
%
%   
%
% Output parameters: Poly_node, a struct that includes the following fields
%   1) success_flag: True if found the unique DSO CRC, False otherwise.
%   2) crc_gen_polys: the final list of CRC candidates; the 
%          list size = 1 (i.e., DSO CRC poly) if success_flag = True, 
%          and list size >1 otherwise;
%   3) stopped_distance: -1 if success_flag = True, and a scalar
%          indicating the stopped distance at which the DSO CRC is
%          still not finalized if success_flag = False;
%   4) crc_distance: a positive integer if undetected errors are
%          found within "d_tilde", and -1 otherise;
%
% Written by Hengjie Yang (hengjie.yang@ucla.edu)   03/30/21
%

tic

Poly_node = {};

if nargin < 7
    base = 16;
end

path = './Simulation_results/';
k = length(numerators);
n = k + 1; 

num_string = '';
for iter = 1:k
    num_string = [num_string, num2str(numerators(iter)), '_'];
end

if mod(K+m, k) ~= 0 % check if (K+m) is divisible by k.
    disp('Invalid parameters of K or m!');
    return
end

fileName = ['error_events_v_',num2str(v), '_num_', num_string, 'den_',...
    num2str(denominator), '_d_tilde_',num2str(d_tilde)];

if ~exist([path, fileName, '.mat'], 'file')
    disp(['Error: the file ',fileName, ' does not exist!']);
    return
end
load([path, fileName, '.mat'], 'error_events', 'error_event_lengths');

error_events = error_events;
error_event_lengths = error_event_lengths;

mu = ceil((v-1)/k); % the length of termination stages
N = (K+m)/k + mu; % the overall trellis length

% fileName = ['weight_spectrum_v_',num2str(v), '_num_', num_string,...
%     'den_', num2str(denominator),'_N_',num2str(N)];
% 
% if ~exist([path, fileName, '.mat'], 'file')
%     disp(['Error: the file ',fileName, ' does not exist!']);
%     return
% end 
% load([path, fileName, '.mat'], 'weight_node');
% weight_spectrum = weight_node.weight_spectrum;
% d_min = find(weight_spectrum > 0);
d_min = 1; % minimum distance


List_size = 2^(m-1);
Candidate_CRCs = dec2bin(0:List_size-1) - '0';
Candidate_CRCs = [ones(List_size,1), Candidate_CRCs, ones(2^(m-1),1)]; % degree order from highest to lowest
Candidate_poly_base = dec2base(bin2dec(num2str(Candidate_CRCs)),base);



Undetected_spectrum = inf(List_size, d_tilde); % each column represents the undected spectrum

if m == 1 % special case
    Candidate_CRCs = [1 1];
    Undetected_spectrum = inf(List_size, 1); % each column represents the undected spectrum
    Candidate_poly_base = dec2base(bin2dec(num2str(Candidate_CRCs)),base);
end

mask = true(size(Candidate_CRCs,1),1);
locations = find(mask == true);

crc_gen_polys = [];
success_flag = false;
stopped_distance = -1;
min_dist = -1;

for dist = d_min:d_tilde
    if ~isempty(error_events{dist})
        weight_vec = zeros(size(locations, 1), 1);

        % Construct single-error events
        parfor ii = 1:size(locations, 1)
            weight_vec(ii) = check_divisible_by_distance(error_events,...
                error_event_lengths, Candidate_CRCs(locations(ii),:), dist, N);
        end

        for ii = 1:size(locations, 1)
            Undetected_spectrum(locations(ii), dist) = weight_vec(ii);
        end

        % Construct double-error events
        if d_tilde >= 2*d_min
            temp = zeros(size(locations, 1), 1);
            parfor ii=1:size(locations,1)
                temp(ii) = check_double_error_divisible_by_distance(error_events, error_event_lengths,...
                    Candidate_CRCs(locations(ii),:), dist, d_min, d_tilde, N);
            end

            for ii = 1:size(locations, 1)
                Undetected_spectrum(locations(ii), dist) = Undetected_spectrum(locations(ii), dist) + temp(ii);
            end
        end

        min_weight = min(weight_vec);
        locations = locations(weight_vec == min_weight);
        disp(['    Current distance: ',num2str(dist),' number of candidates: ',...
            num2str(size(locations,1))]);
        if length(locations) == 1
            crc_gen_polys = Candidate_poly_base(locations(1),:);
            success_flag = true;
            break
        end 
        if dist == d_tilde && length(locations) > 1
            crc_gen_polys = Candidate_poly_base(locations,:);
            stopped_distance = d_tilde;
            disp(['    d_tilde is insufficient to find the DSO CRC... ']);
            disp(['    Stopped distance: ',num2str(stopped_distance),...
                ' # of candidate polynomials: ',num2str(size(crc_gen_polys,1))]);
        end
    end
end





% Identify undetected minimum distance
if success_flag == true
    for dist = d_min:d_tilde
        num = 0;
        if ~isempty(error_events{dist}) 
            num = check_divisible_by_distance(error_events,...
                    error_event_lengths,  Candidate_CRCs(locations(1),:), dist, N);
        end

        if dist>= 2*d_min
            temp = check_double_error_divisible_by_distance(error_events,error_event_lengths,...
                        Candidate_CRCs(locations(1),:), dist, d_min, d_tilde, N);
            num = num + temp;
        end

        if num > 0
            min_dist = dist;
            disp(['    DSO CRC polynomial: ',num2str(crc_gen_polys(1,:))]);
            disp(['    Minimum undetected distance: ',num2str(min_dist)]);
            break
        end      
    end
end


% Save results
Poly_node.success_flag = success_flag;
Poly_node.crc_gen_polys = crc_gen_polys;
Poly_node.stopped_distance = stopped_distance;
Poly_node.crc_distance = min_dist;

fileName = ['Poly_node_v_',num2str(v), '_num_', num_string, 'den_',...
    num2str(denominator), '_K_',num2str(K), '_m_',num2str(m)];
save([path, fileName, '.mat'], 'Poly_node', '-v7.3');


timing = toc;
disp(['Execution time: ',num2str(timing),'s']);

end








function undetected_weight=check_divisible_by_distance(error_event,error_event_length,test_polynomial,dist,max_length)

polynomial = fliplr(test_polynomial); % degree from lowest to highest

undetected_weight = 0;
input = error_event{dist};
input = fliplr(input); % degree from lowest to highest

for i = 1:size(input,1)
    [~,remd] = gfdeconv(input(i,:),polynomial,2);
    if remd == 0
        if error_event_length{dist}(i) <= max_length
            undetected_weight = undetected_weight + max_length - error_event_length{dist}(i) + 1;
        end
    end
end
end



function undetected_error_weight=check_double_error_divisible_by_distance(error_event,error_event_length,test_polynomial,distance,d_min,d_tilde,max_length)

polynomial = fliplr(test_polynomial); %degree from lowest to highest

undetected_error_weight = 0;

for d1 = d_min:d_tilde
    for d2 = d_min:d_tilde
        if d1+d2 == distance && ~isempty(error_event{d1}) && ~isempty(error_event{d2})
            for i = 1:size(error_event{d1},1)
                for j = 1:size(error_event{d2},1)
                    input_1 = fliplr(error_event{d1}(i,:));
                    input_2 = fliplr(error_event{d2}(j,:));%degree from lowest to highest
                    len_1 = error_event_length{d1}(i);
                    len_2 = error_event_length{d2}(j);
                    for g1 = 0:max_length-len_1-len_2
                        double_error_event = [input_2(1:end), zeros(1,g1), input_1(1:end)];
                        [~,remd] = gfdeconv(double_error_event, polynomial, 2);
                        if remd == 0
                            undetected_error_weight = undetected_error_weight+...
                                (max_length-len_1-len_2-g1+1);
                        end
                    end
                end
            end
        end
    end
end

end





