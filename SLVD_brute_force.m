function [check_flag, correct_flag, path_rank] = SLVD_brute_force(rxSig, high_rate_codewords, Psi, input_seq, crc_poly)

% This function is brute-force version of the serial list Viterbi decoder (SLVD).
% The decoder does not require a trellis. It only requires the set of all
% high-rate codewords and emulate the procedure of SLVD.
%
% Input parameters:
%       1) rxSig: the received sequence composed of real numbers
%       2) high_rate_codewords: the set of all high-rate codewords
%       3) Psi: the constrained maximum list rank
%       4) input_seq: the input sequence to the conv. encoder
%       5) crc_poly: the CRC polynomial in binary from lowest to highest
%           degree
%
% Output parameters:
%       1) check_flag: true if SLVD identifies a divisible input sequence
%       within Psi trials.
%       2) correct_flag: true if check_flag is true and the divisible input
%       sequence identified is identical to the transmitted CRC-coded
%       sequence
%       3) path_rank: the terminating list rank of the identified divisible
%       input sequence, a scalar no greater than Psi.
%
% Written by Hengjie Yang (hengjie.yang@ucla.edu)   04/06/21
%

check_flag = -1;
correct_flag = -1;
path_rank = -1;

num_codewords = size(high_rate_codewords, 1);
K = log2(num_codewords);

BPSK_sequences = -2*high_rate_codewords + 1;
Cumulative_metrics = sum((BPSK_sequences - rxSig).^2, 2);

[~, I] = sort(Cumulative_metrics);

cur = 1;
while cur<= Psi
    decoded_input_sequence = I(cur) - 1;
    decoded_input_sequence = dec2bin(decoded_input_sequence, K) - '0';
    [~, remd] = gfdeconv(fliplr(decoded_input_sequence), crc_poly);
    if any(remd) == 0
        check_flag = 1;
        path_rank = cur;
        if all(decoded_input_sequence == input_seq(1:K))
            correct_flag = 1;
        else
            correct_flag = 0;
        end
        return
    end
    cur = cur + 1;
end

check_flag = 0;
path_rank = Psi;

end