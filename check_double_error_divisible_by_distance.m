function undetected_error_weight=check_double_error_divisible_by_distance(error_event,error_event_length,test_polynomial,distance,dfree,d_tilde,max_length)
%
%   Inputs:
%       1) error_event: a d_tilde*1 cell where the i-th cell records all its input
%       2) error_event_length: a d_tilde*1 vector indicating the length of
%               inputs
%       3) test_polynomial: the polynomial to be tested (in binary from highest to lowest degree)
%       4) distance: the specific distance at which we are testing
%       5) d_free: a scalar denoting the free distance of the convolutional code
%       6) d_tilde: a scalar denoting the distance threshold
%       5) max_length: k+m+v
%
%   Outputs:
%       1) undetected_error_weight: a scalar denoting the number of double 
%           undetectable error event.
%
%

%   Copyright 2020 Hengjie Yang



polynomial = fliplr(test_polynomial); %degree from lowest to highest

undetected_error_weight=0;

for d1=dfree:d_tilde
    for d2=dfree:d_tilde
        if d1+d2==distance && ~isempty(error_event{d1}) && ~isempty(error_event{d2})
            for i=1:size(error_event{d1},1)
                for j=1:size(error_event{d2},1)
                    input_1=fliplr(error_event{d1}(i,:));
                    input_2=fliplr(error_event{d2}(j,:));%degree from lowest to highest
                    len_1=error_event_length{d1}(i);
                    len_2=error_event_length{d2}(j);
                    for g1=0:max_length-len_1-len_2
                        double_error_event=[input_2(end-len_2+1:end),zeros(1,g1),input_1(end-len_1+1:end)];
                        [~,remd]=gfdeconv(double_error_event,polynomial,2);
                        if remd==0
                            undetected_error_weight=undetected_error_weight+...
                                (max_length-len_1-len_2-g1+1);
                        end
                    end
                end
            end
        end
    end
end

