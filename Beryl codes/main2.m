clear; clc;
%% Trials with BPSK+AWGN
set(0,'DefaultTextFontName','Times','DefaultTextFontSize',14,...
    'DefaultAxesFontName','Times','DefaultAxesFontSize',14,...
    'DefaultLineLineWidth',1,'DefaultLineMarkerSize',7.75);
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
 
G = [103 0 0 161;
    0 103 0 135;
    0 0 103 155];
H = [1 1 1 0 0 0 1;
    1 0 1 1 1 0 1;
    1 1 0 1 1 0 1;
    1 0 0 0 0 1 1];
n = 4; % rate (n-1)/n
[trellis,mapping] = customtrellis(n,G);
Nblocks = 30; % number of blocks
step = 0.5; % resolution/step size
EbNodB_min = 1; % min snr to test
EbNodB_max = 4; % max snr to test
count = 1; % iteration counter
EbNodB = EbNodB_min:step:EbNodB_max;
errors = zeros(2,size(EbNodB,2)); % stores BER values
Nerrors = 200;
frames = Nblocks/(n-1);
CRC=[1 0 1 1 1; 1 0 0 1 1; 1 1 0 0 1]; %0x17, 0x13, 0x19
for iter = 1:size(EbNodB,2)
    snr = EbNodB(iter);
    disp(['SNR = ',num2str(snr)]);
    count0 = 0; % correct decoded counter for dual list decoder
    count1 = 0; % trial counter for normal Viterbi
    count2 = 0; % error counter for normal Viterbi
    UE=0;  % undetected error counter
    NACK=0; % not pass CRC error counter
    while UE<Nerrors
%         disp(['trials = ',num2str(trial)]);
        % generate message
        msg1 = randi([0,1],1,Nblocks/(n-1));
        msg2 = randi([0,1],1,Nblocks/(n-1));
        msg3 = randi([0,1],1,Nblocks/(n-1));
        crc_msg1 = crc_message(CRC(1,:),msg1);
        crc_msg2 = crc_message(CRC(2,:),msg2);
        crc_msg3 = crc_message(CRC(3,:),msg3);
        crc_msg = [crc_msg1; crc_msg2; crc_msg3];
        crc_msg = reshape(crc_msg,1,[]);
        msg_terminated=find_final_and_terminate(crc_msg,H,mapping);
        % encoder     
        r = convenc(msg_terminated,trellis);
        % BPSK modulation
        txSig = r;
        txSig(txSig == 0) = -1;
        % Send txSig over the AWGN channel
        s = awgn(txSig, snr, 'measured');
        % dual decoder
        dualTrellis = Construct_Dual_Trellis(H);
        L=1e5;
        [check_flag,correct_flag,path_rank,m_cap]=DBS_LVA_Dual_Soft(dualTrellis, n, s, CRC, crc_msg, L);
        % correct decoding
        if correct_flag==1
            count0 = count0+1;
        end
        % undetected error
        if check_flag==1 && ~isequal(m_cap,crc_msg)
            UE = UE+1;
            disp(['Number of Undetected Error = ',num2str(UE)]);
        end
        % NACK
        if check_flag~=1
            NACK=NACK+1;
        end
    end
    errors(1,iter) = UE/(UE+NACK+count0); % FER
%     % Normal Viterbi decoder
%     while count1<(UE+NACK+count0)
% %         disp(['trials = ',num2str(trial)]);
%         % generate message
%         msg = randi([0,1],1,Nblocks);
%         crc_msg = crc_message(CRC,msg);
%         % encoder     
%         r_row = convenc(crc_msg,trellis);
%         r = transpose(reshape(r_row,n,[])); % reshape r_row
%         % BPSK modulation
%         txSig = r;
%         txSig(txSig == 0) = -1;
%         % Send txSig over the AWGN channel
%         s = awgn(txSig, snr, 'measured');
%         % normal Viterbi decoder
%         s_reshape = reshape(s',1,[]);
% %         % Quantize to prepare for soft-decision decoding.
%         demod_s = quantiz(s_reshape,[0.001,.1,.3,.5,.7,.9,.999]);
%         hVitDec = vitdec(demod_s, trellis, 10,'term','soft',3);
% %         temp(count2,2) = sum(abs(crc_msg - hVitDec)) / size(crc_msg,2); % BER
% %         fer(count2,2) = 1-(1-temp(count2,2))^frames; % FER
%         %  UPDATE COUNT ONLY WHEN THERE IS FRAME ERROR
%         if (~isequal(crc_msg,hVitDec))
%             count2 = count2 + 1;
%         end
%         count1=count1+1;
%     end
%     errors(2,iter) = count2/count1; %FER
    UE_NACK(1,iter) = UE/(UE+NACK+count0);
    UE_NACK(2,iter) = NACK/(UE+NACK+count0);
end
 
error_matrix = EbNodB_min:step:EbNodB_max;
figure;
semilogy(error_matrix,errors(1,:),'b');
% set(gca,'xdir','reverse');
grid on
hold on
% semilogy(error_matrix,errors(2,:),'r');
title('FER vs. SNR for r=3/4 Terminated Convolutional Code');
xlabel('SNR');
ylabel('Frame Error Rate (FER)');
% legend('Dual List Decoding','Normal Viterbi Decoding');
 
figure;
semilogy(error_matrix,UE_NACK(1,:),'b');
set(gca,'xdir','reverse');
grid on
hold on
semilogy(error_matrix,UE_NACK(2,:),'r');
title('Error Probability vs. SNR for Undetected Error and NACK');
xlabel('SNR');
ylabel('Probability');
legend('Undetected Error','NACK');
