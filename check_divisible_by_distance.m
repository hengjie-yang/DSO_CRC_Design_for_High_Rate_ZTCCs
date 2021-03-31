function undetected_weight=check_divisible_by_distance(error_event,error_event_length,test_polynomial,d,max_length)

%
%   Inputs:
%       1) error_event: a d_tilde*1 cell where the i-th cell records all its input
%       2) test_polynomial: the polynomial to be tested (in binary from highest to lowest degree)
%       3) d: the specific distance at which we are testing
%       4) max_length: overall trellis length
%
%   Outputs:
%       1) undetected_weight: a scalar indicating the number of
%       undetectable single error events of distance 'd'
%

%   Copyright 2020 Hengjie Yang

polynomial = fliplr(test_polynomial); % degree from lowest to highest

undetected_weight = 0;
input = error_event{d};
input = fliplr(input); % degree from lowest to highest

for i = 1:size(input,1)
    [~,remd] = gfdeconv(input(i,:),polynomial,2);
    if remd == 0 
        if error_event_length{d}(i) <= max_length
            undetected_weight = undetected_weight + max_length - error_event_length{d}(i) + 1;
        end
    end
end
