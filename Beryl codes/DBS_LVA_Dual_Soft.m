function [check_flag,correct_flag,path_rank,dec]=DBS_LVA_Dual_Soft(trellis, n, rxSig, polynomial, input, L)
%Input: 1)trellis: the CC trellis; 
%       2)rxSig: the received signal
%       3)polynomial: the CRC code, in binary; 
%       4)input: the original input sequence (9 bits)
%       5)L: maximum list size

%Output: 1)check_flag: whether LVA finds a codeword that checks CRC
%        2)correct_flag: whether the codeword found is equal to the input
%        3)path_rank: the codeword rank, max: Trial_threshold
%        4)dec: the decoded input

%Notes: 
%       1)In simulation, the input has to be global variable to obtain correct_flag
%       2)The program works for rate-1/omega convolutional code.

global NumInput
NumInput=trellis.numInputSymbols;
NumOutBit=trellis.numOutputSymbols;
NumState=trellis.numStates;
% v=log2(NumState);
length=size(rxSig,2)/NumOutBit;%message length
k=size(input,2)/(n-1);

global Column
Column=cell(length,1);


for iteration=0:length-1
    Column{iteration+1}=cell(NumState,1);%initialization all states for next depth
    received_point=rxSig(NumOutBit*iteration+1:NumOutBit*iteration+NumOutBit);
    
    %update next_state
    if iteration==0
        for i=1:NumInput
            temp=trellis.nextStates(1,1,i);
            if temp ~= Inf
                next_state=temp;
                output=i-1;
                output=dec2bin(oct2dec(output), NumOutBit) - '0'; 
                output_point=Get_Point(output);
                diff_metric=sum((received_point-output_point).^2);
%                 diff_metric=sum(received_point~=output);%Hamming distance
            end
            if isempty(Column{iteration+1}{next_state})
                Column{iteration+1}{next_state}=cell(2,1);%minimum_metric, father state
                Column{iteration+1}{next_state}{1}=[diff_metric,iteration+1];
                Column{iteration+1}{next_state}{2}=[Inf, -1];%[difference, state]
            end
        end
%     elseif iteration<length-v %message part
    else % iteration<length %message part
        for i=1:NumState
            if ~isempty(Column{iteration}{i})
                for j=1:NumInput
                    temp=trellis.nextStates(i,mod(iteration,n)+1,j);
                    if temp ~= Inf
                        next_state=temp;
                        output=j-1; 
                        output=dec2bin(oct2dec(output),NumOutBit) - '0';
                        output_point=Get_Point(output);
                        diff_metric=sum((received_point-output_point).^2)+Column{iteration}{i}{1}(1);
%                         diff_metric=sum(received_point~=output)+Column{iteration}{i}{1}(1);
                        temp1=Column{iteration+1}{next_state};
                        if isempty(temp1)
                            Column{iteration+1}{next_state}=cell(2,1);%minimum_metric, father state
                            Column{iteration+1}{next_state}{1}=[diff_metric,i];%optimal path metric, father state
                            Column{iteration+1}{next_state}{2}=[Inf, -1];%differential metric, father state
                        elseif Column{iteration+1}{next_state}{1}(1)>diff_metric%already exist a better metric
                            diff=Column{iteration+1}{next_state}{1}(1)-diff_metric;
                            state=Column{iteration+1}{next_state}{1}(2);
                            Column{iteration+1}{next_state}{1}=[diff_metric,i];
                            Column{iteration+1}{next_state}{2}=[diff,state];
                        else %two paths going into next state, but they have the same metric
                            %%only one path going into next state, and this path is suboptimal
%                             diff=diff_metric-Column{iteration+1}{next_state}{1}(1);
%                             Column{iteration+1}{next_state}{2}=[diff,i];
                            Column{iteration+1}{next_state}{2}=[diff_metric,i];
                        end
                    end
                end
            end
        end
    end
end


decoded=[];
global Detour_Tree
global node_cnt
global Detour_Array
global RowToID_Table
% global Max
% Max=L+1;
Detour_Tree=[];
node_cnt=1;
loop_cnt=1;
Node=struct('depth',length+5,'father',-1);
Detour_Tree=[Detour_Tree;Node];
Detour_Array=[];
RowToID_Table=[];
q=PriorityQueue([1, -1]);
diff_metric=0;
counter=0;
path_rank=L;
check_flag=-1;
correct_flag=-1;
% distance_diff=-1;
% input_metric=Evaluate_metric_Hamming(trellis, rxSig, input);
dec=-1;
% while abs(diff_metric)<Metric_thrld
all_codewords=[];
while counter<L
    counter=counter+1;
%     disp(['Current trial: ',num2str(counter)]);
    
    node=q.remove();
    id=node(1);
    [codeword, next_state_ids, GMD]=trace_back_search(trellis,length,id);
    c_cap=reshape(codeword,n,[])';
    codeword=find_original_message_v3(c_cap,k);
    all_codewords=[all_codewords; codeword];
%     disp(codeword);
    decoded=[decoded;codeword];
    if next_state_ids(1)~=-1
        for i=1:size(next_state_ids,2)
            q.insert([next_state_ids(i), GMD]);
        end
    end
%     codeword_metric=Evaluate_metric_Hamming(trellis,rxSig, codeword);
%     diff_metric=codeword_metric-input_metric;%Compute the metric difference
    
    %CRC check
%     dataSoft=fliplr(codeword);%degree from low to high
    dataSoft=codeword;
    dataSoft=reshape(dataSoft,(n-1),[]);
    remd=[];
    for row=1:n-1
        [~,temp]=gfdeconv(dataSoft(row,:),polynomial(row,:),2);
        remd=[remd temp];
    end
%     [~,remd]=gfdeconv(dataSoft,polynomial,2);
    if any(remd)==0
        check_flag=1;
        path_rank=counter;
%         distance_diff=sum(codeword~=input);
        dec=codeword;
        if all(codeword==input)
            correct_flag=1;
            break
        else
            correct_flag=0;
        end
% %         disp(all_codewords);
        return %return the whole function
    end
end
% check_flag=0;
% disp(size(all_codewords,1)-1);
% distance_diff=sum(input~=decoded(end,:));
    
    
        

%trace back to find all the L-best paths
function [codeword, next_state_ids, GMD_value]=trace_back_search(trellis, length, state_id)
%Input parameters:
%   1)length: the codeword length
%   2)detour_depth: the depth where we switch the trace trajectory
%   3)state_id: the leaf ID of the path, which is a row vector
%Output parameters:
%   1)codeword: the decoded message
%   2)next_state_ids: the next states with current minimum GMD
%   3)GMD_value: current minimum GMD 
global Column
global NumInput
global Detour_Array
global Detour_Tree
global RowToID_Table

n=size(trellis.nextStates,2);
detour_paths=DFS(state_id);
codeword=[];
GMD=0;
LMD=zeros(1,length);
current_state=1;
for iteration=length:-1:1
%     temp=Column{iteration}{current_state}{2}(1);
%     if temp ~= Inf
    LMD(iteration)=Column{iteration}{current_state}{2}(1)+GMD;
%     end
    if any(detour_paths==iteration)%iteration is in the list, switch direction
        pre_state=Column{iteration}{current_state}{2}(2);
        GMD=GMD+Column{iteration}{current_state}{2}(1);
    else
        pre_state=Column{iteration}{current_state}{1}(2);
    end
    for i=1:NumInput
        if trellis.nextStates(pre_state,mod((iteration-1),n)+1,i)==current_state
            codeword=[codeword,i-1];
            current_state=pre_state;
            break
        end
    end
end
threshold=min(detour_paths);
LMD(threshold:end)=Inf;
Detour_Array=[Detour_Array;LMD];
RowToID_Table(end+1)=state_id;
codeword=codeword(end:-1:1);
GMD_value=min(min(Detour_Array));
if GMD_value==Inf %Sometimes this may happen
    next_state_ids=-1;
    return
end
[rows,next_depths]=find(Detour_Array==GMD_value);
Detour_Array(rows,next_depths)=Inf; 
%repace with Inf so that we can find next minimum value

%Map rows to leaf IDs on the detour tree
rows=rows(:);
next_depths=next_depths(:);
father_ids=RowToID_Table(rows);

next_state_ids=Add_nodes(father_ids,next_depths);
% disp(['Current_min_metric: ',num2str(min_metric)]);




function detour_paths=DFS(state_id)
%state_id: the leaf_state
%the function is to find the path from leaf_state to the root
global Detour_Tree
% global Max
detour_paths=[];
path_rank=1;
cur=state_id;
if cur==1
    return
end
while cur~=1
    detour_paths=[detour_paths,Detour_Tree(cur).depth];
    cur=Detour_Tree(cur).father;
end
% for i=1:Max
%     depth=Detour_Tree(i).depth;
%     if depth==detour_paths(1)
%         path_rank=i;
%         break
%     end
% end
detour_paths=detour_paths(end:-1:1);

   
function next_state_ids=Add_nodes(path_ranks,depths)
%path_rank: current path
%depths: the newly found states on the tree
%loop_cnt: current loop number
global Detour_Tree
global node_cnt
next_state_ids=[];
for i=1:size(depths,1)
    node_cnt=node_cnt+1;
    node=struct('depth',depths(i),'father',path_ranks(i));
    Detour_Tree=[Detour_Tree;node];
    next_state_ids=[next_state_ids,node_cnt];
end

function point=Get_Point(n)
if n==0
    point=-1;
elseif n==1
    point=1;
end

