function [K_num]=KJZ(num,K,design_num)

K_num=zeros(1,design_num);
for i=1:1:design_num
    K_num(i)=mod(num,K);
    num=(num-K_num(i))/K;
end

