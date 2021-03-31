function m_cap = find_original_message_v3(c_cap,k) 
% for v=6 optimal rate 3/4 H(D)
% c_cap: num total blocks*n
% k is number of information blocks
m_cap=[];
n=size(c_cap,2);
for i=1:k
    temp=c_cap(i,1:n-1);
    m_cap=[m_cap;temp];
end
m_cap=reshape(m_cap',[],1)';