function [tau,omega,mu,sigma,MU,SIGMA,z,v,log_likelihood]=MLE(sample,K,L,design_num,context_num)
%initialize tau
tau=zeros(1,K);
for k=1:1:K
    tau(k)=k;
end
tau=tau/sum(tau);
%initialize omega
omega=zeros(1,L);
for l=1:1:L
    omega(l)=l;
end
omega=omega/sum(omega);
%initialize MU
MU=zeros(K,L);
for k=1:1:K
    for l=1:1:L
        MU(k,l)=100*((k-1)*l+l)/K/L;
    end
end
%initialize SIGMA
SIGMA=ones(K,L)*3;
%sample_sigma
sample_sigma=zeros(design_num,context_num);
for i=1:1:design_num
    for j=1:1:context_num
        sample_sigma(i,j)=10;
    end
end
%compute mu sigma C
mu=zeros(design_num,context_num,K,L);
sigma=zeros(design_num,context_num,K,L);
C=zeros(design_num,context_num,K,L);
for i=1:1:design_num
    for j=1:1:context_num
        for k=1:1:K
            for l=1:1:L
                sigma(i,j,k,l)=1/sqrt(length(sample{i,j})/(sample_sigma(i,j)^2)+1/(SIGMA(k,l)^2));
                mu(i,j,k,l)=(sigma(i,j,k,l)^2)*(sum(sample{i,j})/(sample_sigma(i,j)^2)+MU(k,l)/(SIGMA(k,l)^2));
                C(i,j,k,l)=log((1/(2*pi*(sample_sigma(i,j)^2)))^(length(sample{i,j})/2)*sigma(i,j,k,l)/SIGMA(k,l)*exp(((mu(i,j,k,l)/sigma(i,j,k,l))^2-sum(sample{i,j}.^2)/(sample_sigma(i,j)^2)-(MU(k,l)/SIGMA(k,l))^2)/2));
            end
        end
    end
end
%compute z
temp=zeros(K^design_num,L^context_num);
for num1=0:1:(K^design_num-1)
    temp_designcluster=KJZ(num1,K,design_num)+ones(1,design_num);
    for num2=0:1:(L^context_num-1)
        temp_contextcluster=KJZ(num2,L,context_num)+ones(1,context_num);
        for i=1:1:design_num
            temp(num1+1,num2+1)=temp(num1+1,num2+1)+log(tau(temp_designcluster(i)));
        end
        for j=1:1:context_num
            temp(num1+1,num2+1)=temp(num1+1,num2+1)+log(omega(temp_contextcluster(j)));
        end
        for i=1:1:design_num
            for j=1:1:context_num
                temp(num1+1,num2+1)=temp(num1+1,num2+1)+C(i,j,temp_designcluster(i),temp_contextcluster(j));
            end
        end
    end
end
temp=temp-ones(K^design_num,L^context_num)*max(max(temp));
z=zeros(design_num,K);
for i=1:1:design_num
    for k=1:1:K
        for num1=0:1:(K^design_num-1)
            temp_designcluster=KJZ(num1,K,design_num)+ones(1,design_num);
            if temp_designcluster(i)==k
               z(i,k)=z(i,k)+sum(exp(temp(num1+1,:)))/sum(sum(exp(temp)));
            end 
        end
    end
end
%compute v
v=zeros(context_num,L);
for j=1:1:context_num
    for l=1:1:L
        for num2=0:1:(L^context_num-1)
            temp_contextcluster=KJZ(num2,L,context_num)+ones(1,context_num);
            if temp_contextcluster(j)==l
               v(j,l)=v(j,l)+sum(exp(temp(:,num2+1)))/sum(sum(exp(temp)));
            end 
        end
    end
end
%update
for s=1:1:1
    %update tau
    for k=1:1:K
        tau(k)=sum(z(:,k))/design_num;
    end
    %update omega
    for l=1:1:L
        omega(l)=sum(v(:,l))/context_num;
    end
    %update MU SIGMA
    for k=1:1:K
        for l=1:1:L
            temp1=0;
            temp2=0;
            temp3=0;
            for i=1:1:design_num
                for j=1:1:context_num
                    temp1=temp1+z(i,k)*v(j,l);
                    temp2=temp2+z(i,k)*v(j,l)*mu(i,j,k,l);
                end
            end  
            MU(k,l)=temp2/temp1;
            for i=1:1:design_num
                for j=1:1:context_num
                    temp3=temp3+z(i,k)*v(j,l)*(sigma(i,j,k,l)^2+(mu(i,j,k,l)-MU(k,l))^2);
                end
            end
            SIGMA(k,l)=sqrt(temp3/temp1);
        end
    end
    %update mu sigma C
    for i=1:1:design_num
        for j=1:1:context_num
            for k=1:1:K
                for l=1:1:L
                    sigma(i,j,k,l)=1/sqrt(length(sample{i,j})/(sample_sigma(i,j)^2)+1/(SIGMA(k,l)^2));
                    mu(i,j,k,l)=(sigma(i,j,k,l)^2)*(sum(sample{i,j})/(sample_sigma(i,j)^2)+MU(k,l)/(SIGMA(k,l)^2));
                    C(i,j,k,l)=log((1/(2*pi*(sample_sigma(i,j)^2)))^(length(sample{i,j})/2)*sigma(i,j,k,l)/SIGMA(k,l)*exp(((mu(i,j,k,l)/sigma(i,j,k,l))^2-sum(sample{i,j}.^2)/(sample_sigma(i,j)^2)-(MU(k,l)/SIGMA(k,l))^2)/2));
                end
            end
        end
    end
    %update z
    temp=zeros(K^design_num,L^context_num);
    for num1=0:1:(K^design_num-1)
        temp_designcluster=KJZ(num1,K,design_num)+ones(1,design_num);
        for num2=0:1:(L^context_num-1)
            temp_contextcluster=KJZ(num2,L,context_num)+ones(1,context_num);
            for i=1:1:design_num
                temp(num1+1,num2+1)=temp(num1+1,num2+1)+log(tau(temp_designcluster(i)));
            end
            for j=1:1:context_num
                temp(num1+1,num2+1)=temp(num1+1,num2+1)+log(omega(temp_contextcluster(j)));
            end
            for i=1:1:design_num
                for j=1:1:context_num
                    temp(num1+1,num2+1)=temp(num1+1,num2+1)+C(i,j,temp_designcluster(i),temp_contextcluster(j));
                end
            end
        end
    end
    temp=temp-ones(K^design_num,L^context_num)*max(max(temp));
    z=zeros(design_num,K);
    for i=1:1:design_num
        for k=1:1:K
            for num1=0:1:(K^design_num-1)
                temp_designcluster=KJZ(num1,K,design_num)+ones(1,design_num);
                if temp_designcluster(i)==k
                   z(i,k)=z(i,k)+sum(exp(temp(num1+1,:)))/sum(sum(exp(temp)));
                end 
            end
        end
    end
    %update v
    v=zeros(context_num,L);
    for j=1:1:context_num
        for l=1:1:L
            for num2=0:1:(L^context_num-1)
                temp_contextcluster=KJZ(num2,L,context_num)+ones(1,context_num);
                if temp_contextcluster(j)==l
                   v(j,l)=v(j,l)+sum(exp(temp(:,num2+1)))/sum(sum(exp(temp)));
                end 
            end
        end
    end
end
for num1=0:1:(K^design_num-1)
    temp_designcluster=KJZ(num1,K,design_num)+ones(1,design_num);
    for num2=0:1:(L^context_num-1)
        temp_contextcluster=KJZ(num2,L,context_num)+ones(1,context_num);
        for i=1:1:design_num
            temp(num1+1,num2+1)=temp(num1+1,num2+1)+log(tau(temp_designcluster(i)));
        end
        for j=1:1:context_num
            temp(num1+1,num2+1)=temp(num1+1,num2+1)+log(omega(temp_contextcluster(j)));
        end
        for i=1:1:design_num
            for j=1:1:context_num
                temp(num1+1,num2+1)=temp(num1+1,num2+1)+C(i,j,temp_designcluster(i),temp_contextcluster(j));
            end
        end
    end
end
log_likelihood=max(max(temp));
temp=temp-ones(K^design_num,L^context_num)*max(max(temp));
log_likelihood=log_likelihood+log(sum(sum(exp(temp))));


            