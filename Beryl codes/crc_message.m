function crc_msg = crc_message(CRC,msg)
msg_shifted = [msg zeros(1,(size(CRC,2)-1))]; %shift msg by 6 bits (degree of CRC)
r=zeros(1,(size(CRC,2)-1));
%divide each shifted frame of msg by CRC
[~,remainder] = gfdeconv(fliplr(msg_shifted),fliplr(CRC));
% [~,remainder] = gfdeconv(msg,CRC);
r(1,(end-size(remainder,2)+1):end)=fliplr(remainder);
crc_msg = [msg'; r']'; %concatenate CRC to each frame of the message
