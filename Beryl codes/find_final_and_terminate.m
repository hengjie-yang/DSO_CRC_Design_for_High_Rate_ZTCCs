function msg_terminated = find_final_and_terminate(crc_msg,H,mapping)
% find final state and determine termination bits
msg_reshape=reshape(crc_msg',3,[]);
u1=msg_reshape(1,:);
u2=msg_reshape(2,:);
u3=msg_reshape(3,:);
p1=gfconv(fliplr(u1),fliplr(H(1,:)));
p2=gfconv(fliplr(u2),fliplr(H(2,:)));
p3=gfconv(fliplr(u3),fliplr(H(3,:)));
p_all=gfadd(p3,gfadd(p1,p2));
[~,rmd]=gfdeconv(p_all,fliplr(H(4,:)));
rmdsize=size(rmd,2);

if rmdsize<6
    rmd=[rmd zeros(1,(6-rmdsize))];
end
rmd10=bi2de(fliplr(rmd));

[bits1,bits2]=find(mapping==rmd10);
bits12=de2bi(bits1-1,3);
bits22=de2bi(bits2-1,3);

termination_bits=fliplr([bits12;bits22])';
msg_terminated=[msg_reshape termination_bits];
% m1=msg_terminated(1,:);
% m2=msg_terminated(2,:);
% m3=msg_terminated(3,:);
% curstate=0;
% for i=1:size(m1,2)
%     input=[m1(1,i) m2(1,i) m3(1,i)];
%     input = bi2de(fliplr(input));
%     nextstate=trellis.nextStates(curstate+1,input+1);
%     curstate=nextstate;
% end

msg_terminated=reshape(msg_terminated,1,[]);