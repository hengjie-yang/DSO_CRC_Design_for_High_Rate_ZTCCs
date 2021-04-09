function [check_flag, correct_flag, path_rank, dec] = SLVD_dual_trellis(dual_trellis, rxSig, crc_poly, tx_seq, Psi)

%
% This function is to perform serial list Viterbi decoding (SLVD) on the dual
% trellis of a (n, n-1, v) convolutional encoder. The detailed dual trellis
% construction can be found in [1].
%
% References:
%   [1] T. Yamada et al., "A new maximum likelihood decoding of high rate
%       convolutional codes using a trellis", Elec. Commun. in Japan, 1983.
%
% Input parameters:
%   1) dual_trellis, a struct denoting the overall dual trellis
%   2) rxSig: a length n*h vector denoting the received real signals, where
%           h is the number of input-output streams.
%   3) crc_poly: CRC polynonmial in binary, degree from lowest to highest
%   4) tx_seq: a length (n-1)*h binary vector denoting the input sequence
%           to the (n, n-1, v) convolutional encoder. This is only to used
%           to determine "correct_flag".
%   5) Psi: a scalar denoting the constrained maximum list size.
%
% Output parameters:
%   1) check_flag: true if the decoded input sequence is divisible by
%           "crc_poly".
%   2) correct_flag: true if the decoded input sequence is indentical to
%           "tx_seq".
%   3) path_rank: a scalar denoting the terminating list rank at which
%           divisibility is first met.
%   4) dec: the identified decoded input sequence.
%
% Remarks:
%   1) The decoder only works for (n, n-1, v) conv. encoder.
%
% Written by Hengjie Yang (hengjie.yang@ucla.edu)   04/08/21.
%



















