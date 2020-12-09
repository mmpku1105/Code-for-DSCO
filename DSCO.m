function [CS]=DSCO(sim_model,design_num,context_num,K,L,T,true_best)  %%sim_model is any simulation model

%initialize sample
sample=cell(design_num,context_num);
initial_k=5;
for i=1:1:design_num
    for j=1:1:context_num
        for k=1:1:initial_k
            sample{i,j}=[sample{i,j},sim_model(i,j)];
        end
    end
end

[~,~,mu,sigma,MU,SIGMA,z,v,~]=MLE(sample,K,L,design_num,context_num);

%performance estimation
estimate_y=zeros(design_num,context_num);
for i=1:1:design_num
    for j=1:1:context_num
        [~,temp_k]=max(z(i,:));
        [~,temp_l]=max(v(j,:));
        estimate_y(i,j)=mu(i,j,temp_k,temp_l);
    end
end

%compute VFA
VFA=0;
best_design=zeros(1,context_num);
for j=1:1:context_num
    temp=0;
    [~,best_design(j)]=max(estimate_y(:,j));
    for num1=0:1:(K^design_num-1)
        temp_designcluster=KJZ(num1,K,design_num)+ones(1,design_num);
        for l=1:1:L
            temp1=v(j,l);
            temp2=0;
            for i=1:1:design_num
                temp1=temp1*z(i,temp_designcluster(i));
                if (i~=best_design(j))
                    temp3=(mu(best_design(j),j,temp_designcluster(best_design(j)),l)-mu(i,j,temp_designcluster(i),l))^2/(sigma(best_design(j),j,temp_designcluster(best_design(j)),l)^2+sigma(i,j,temp_designcluster(i),l)^2);
                    if temp2==0
                        temp2=temp3;
                    else
                        temp2=min(temp2,temp3);
                    end
                end
            end
            temp1=temp1*temp2;
            temp=temp+temp1;
        end
    end
    if j==1
        VFA=temp;
    else
        VFA=min(VFA,temp);
    end
end

%% sequential sampling
t=1;
while(t<=T)
    VFA_onestep=zeros(design_num,context_num);
    for r=1:1:design_num
        for q=1:1:context_num
            [~,temp_k]=max(z(r,:));
            [~,temp_l]=max(v(q,:));
            z_onestep=z;
            for k=1:1:K
                z_onestep(r,k)=log(z(r,k)+0.00001)-(mu(r,q,k,temp_l)-MU(k,temp_l))^2*sigma(r,q,k,temp_l)^4/(2*sample_sigma(r,q)^2*SIGMA(k,temp_l)^4);
            end
            z_onestep(r,:)=z_onestep(r,:)-ones(1,K)*max(z_onestep(r,:));
            z_onestep(r,:)=exp(z_onestep(r,:))/sum(exp(z_onestep(r,:)));
            v_onestep=v;
            for l=1:1:L
                v_onestep(q,l)=log(v(q,l)+0.00001)-(mu(r,q,temp_k,l)-MU(temp_k,l))^2*sigma(r,q,temp_k,l)^4/(2*sample_sigma(r,q)^2*SIGMA(temp_k,l)^4);
            end
            v_onestep(q,:)=v_onestep(q,:)-ones(1,L)*max(v_onestep(q,:));
            v_onestep(q,:)=exp(v_onestep(q,:))/sum(exp(v_onestep(q,:)));
            SIGMA_onestep=SIGMA;
            for k=1:1:K
                for l=1:1:L
                    temp1=0;
                    temp2=0;
                    for i=1:1:design_num
                        for j=1:1:context_num
                            temp1=temp1+z(i,k)*v(j,l);
                            if (i==r)&&(j==q)
                                temp2=temp2+z(i,k)*v(j,l)*(1/((length(sample{i,j})+1)/(sample_sigma(i,j)^2)+1/SIGMA(k,l)^2)+(mu(i,j,k,l)-MU(k,l))^2);
                            else
                                temp2=temp2+z(i,k)*v(j,l)*(1/(length(sample{i,j})/(sample_sigma(i,j)^2)+1/SIGMA(k,l)^2)+(mu(i,j,k,l)-MU(k,l))^2);
                            end
                        end
                    end
                    SIGMA_onestep(k,l)=sqrt(temp2/temp1);
                end
            end
            for j=1:1:context_num
                temp=0;
                for num1=0:1:(K^design_num-1)
                    temp_designcluster=KJZ(num1,K,design_num)+ones(1,design_num);
                    for l=1:1:L
                        sigma_onestep=sigma;
                        for i=1:1:design_num
                            if (i==r)&&(j==q)
                                sigma_onestep(i,j,temp_designcluster(i),l)=sqrt(1/((length(sample{i,j})+1)/(sample_sigma(i,j)^2)+1/SIGMA_onestep(temp_designcluster(i),l)^2));
                            elseif (temp_designcluster(i)==temp_k)&&(l==temp_l)
                                sigma_onestep(i,j,temp_designcluster(i),l)=sqrt(1/(length(sample{i,j})/(sample_sigma(i,j)^2)+1/SIGMA_onestep(temp_designcluster(i),l)^2));
                            end
                        end
                        temp1=v_onestep(j,l);
                        temp2=0;
                        for i=1:1:design_num
                            temp1=temp1*z_onestep(i,temp_designcluster(i));
                            if (i~=best_design(j))
                                temp3=(mu(best_design(j),j,temp_designcluster(best_design(j)),l)-mu(i,j,temp_designcluster(i),l))^2/(sigma_onestep(best_design(j),j,temp_designcluster(best_design(j)),l)^2+sigma_onestep(i,j,temp_designcluster(i),l)^2);
                                if temp2==0
                                    temp2=temp3;
                                else
                                    temp2=min(temp2,temp3);
                                end
                            end
                        end
                        temp1=temp1*temp2;
                        temp=temp+temp1;
                    end
                end
                if j==1
                    VFA_onestep(r,q)=temp;
                else
                    VFA_onestep(r,q)=min(VFA_onestep(r,q),temp);
                end
            end
        end
    end
    if max(max(VFA_onestep))>VFA
        [r_set,q_set]=find(VFA_onestep==max(max(VFA_onestep)));
        r=r_set(end);
        q=q_set(end);
        sample{r,q}=[sample{r,q},sim_model(r,q)];
        t=t+1;
        %update
        [~,temp_k]=max(z(r,:));
        [~,temp_l]=max(v(q,:));
        for k=1:1:K
            z(r,k)=log(z(r,k)+0.00001)-(mu(r,q,k,temp_l)-MU(k,temp_l))^2*sigma(r,q,k,temp_l)^4/(2*sample_sigma(r,q)^2*SIGMA(k,temp_l)^4);
        end
        z(r,:)=z(r,:)-ones(1,K)*max(z(r,:));
        z(r,:)=exp(z(r,:))/sum(exp(z(r,:)));
        for l=1:1:L
            v(q,l)=log(v(q,l)+0.00001)-(mu(r,q,temp_k,l)-MU(temp_k,l))^2*sigma(r,q,temp_k,l)^4/(2*sample_sigma(r,q)^2*SIGMA(temp_k,l)^4);
        end
        v(q,:)=v(q,:)-ones(1,L)*max(v(q,:));
        v(q,:)=exp(v(q,:))/sum(exp(v(q,:)));
        for k=1:1:K
            for l=1:1:L
                sigma(r,q,k,l)=1/sqrt(length(sample{r,q})/(sample_sigma(r,q)^2)+1/(SIGMA(k,l)^2));
                mu(r,q,k,l)=(sigma(r,q,k,l)^2)*(sum(sample{r,q})/(sample_sigma(r,q)^2)+MU(k,l)/(SIGMA(k,l)^2));
            end
        end
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
        for k=1:1:K
            for l=1:1:L
                sigma(r,q,k,l)=1/sqrt(length(sample{r,q})/(sample_sigma(r,q)^2)+1/(SIGMA(k,l)^2));
                mu(r,q,k,l)=(sigma(r,q,k,l)^2)*(sum(sample{r,q})/(sample_sigma(r,q)^2)+MU(k,l)/(SIGMA(k,l)^2));
            end
        end
    else
        W=0;
        while((W<VFA)&&(t<=T))
            W_onestep=zeros(design_num,context_num);
            for r=1:1:design_num
                for q=1:1:context_num
                    SIGMA_onestep=SIGMA;
                    for k=1:1:K
                        for l=1:1:L
                            temp1=0;
                            temp2=0;
                            for i=1:1:design_num
                                for j=1:1:context_num
                                    temp1=temp1+z(i,k)*v(j,l);
                                    if (i==r)&&(j==q)
                                        temp2=temp2+z(i,k)*v(j,l)*(1/((length(sample{i,j})+1)/(sample_sigma(i,j)^2)+1/SIGMA(k,l)^2)+(mu(i,j,k,l)-MU(k,l))^2);
                                    else
                                        temp2=temp2+z(i,k)*v(j,l)*(1/(length(sample{i,j})/(sample_sigma(i,j)^2)+1/SIGMA(k,l)^2)+(mu(i,j,k,l)-MU(k,l))^2);
                                    end
                                end
                            end
                            SIGMA_onestep(k,l)=sqrt(temp2/temp1);
                        end
                    end
                    for j=1:1:context_num
                        temp1=0;
                        for num1=0:1:(K^design_num-1)
                            temp_designcluster=KJZ(num1,K,design_num)+ones(1,design_num);
                            for l=1:1:L
                                sigma_onestep=sigma;
                                for i=1:1:design_num
                                    if (i==r)&&(j==q)
                                        sigma_onestep(i,j,temp_designcluster(i),l)=sqrt(1/((length(sample{i,j})+1)/(sample_sigma(i,j)^2)+1/SIGMA_onestep(temp_designcluster(i),l)^2));
                                    elseif (temp_designcluster(i)==temp_k)&&(l==temp_l)
                                        sigma_onestep(i,j,temp_designcluster(i),l)=sqrt(1/(length(sample{i,j})/(sample_sigma(i,j)^2)+1/SIGMA_onestep(temp_designcluster(i),l)^2));
                                    end
                                end
                                temp2=0;
                                for i=1:1:design_num
                                    if (i~=best_design(j))
                                        temp3=(mu(best_design(j),j,temp_designcluster(best_design(j)),l)-mu(i,j,temp_designcluster(i),l))^2/(sigma_onestep(best_design(j),j,temp_designcluster(best_design(j)),l)^2+sigma_onestep(i,j,temp_designcluster(i),l)^2);
                                        if temp2==0
                                            temp2=temp3;
                                        else
                                            temp2=min(temp2,temp3);
                                        end
                                    end
                                end
                                if temp1==0
                                    temp1=temp2;
                                else
                                    temp1=min(temp1,temp2);
                                end
                            end
                        end
                        if j==1
                            W_onestep(r,q)=temp1;
                        else
                            W_onestep(r,q)=min(W_onestep(r,q),temp1);
                        end
                    end
                end
            end
            W=max(max(W_onestep));
            [r_set,q_set]=find(W_onestep==max(max(W_onestep)));
            r=r_set(end);
            q=q_set(end);
            sample{r,q}=[sample{r,q},sim_model(r,q)];
            t=t+1;
            %update
            [~,temp_k]=max(z(r,:));
            [~,temp_l]=max(v(q,:));
            for k=1:1:K
                z(r,k)=log(z(r,k)+0.00001)-(mu(r,q,k,temp_l)-MU(k,temp_l))^2*sigma(r,q,k,temp_l)^4/(2*sample_sigma(r,q)^2*SIGMA(k,temp_l)^4);
            end
            z(r,:)=z(r,:)-ones(1,K)*max(z(r,:));
            z(r,:)=exp(z(r,:))/sum(exp(z(r,:)));
            for l=1:1:L
                v(q,l)=log(v(q,l)+0.00001)-(mu(r,q,temp_k,l)-MU(temp_k,l))^2*sigma(r,q,temp_k,l)^4/(2*sample_sigma(r,q)^2*SIGMA(temp_k,l)^4);
            end
            v(q,:)=v(q,:)-ones(1,L)*max(v(q,:));
            v(q,:)=exp(v(q,:))/sum(exp(v(q,:)));
            for k=1:1:K
                for l=1:1:L
                    sigma(r,q,k,l)=1/sqrt(length(sample{r,q})/(sample_sigma(r,q)^2)+1/(SIGMA(k,l)^2));
                    mu(r,q,k,l)=(sigma(r,q,k,l)^2)*(sum(sample{r,q})/(sample_sigma(r,q)^2)+MU(k,l)/(SIGMA(k,l)^2));
                end
            end
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
            for k=1:1:K
                for l=1:1:L
                    sigma(r,q,k,l)=1/sqrt(length(sample{r,q})/(sample_sigma(r,q)^2)+1/(SIGMA(k,l)^2));
                    mu(r,q,k,l)=(sigma(r,q,k,l)^2)*(sum(sample{r,q})/(sample_sigma(r,q)^2)+MU(k,l)/(SIGMA(k,l)^2));
                end
            end
        end
    end
    %performance estimation
    estimate_y=zeros(design_num,context_num);
    for i=1:1:design_num
        for j=1:1:context_num
            [~,temp_k]=max(z(i,:));
            [~,temp_l]=max(v(j,:));
            estimate_y(i,j)=mu(i,j,temp_k,temp_l);
        end
    end
    %compute VFA
    VFA=0;
    best_design=zeros(1,context_num);
    for j=1:1:context_num
        temp=0;
        [~,best_design(j)]=max(estimate_y(:,j));
        for num1=0:1:(K^design_num-1)
            temp_designcluster=KJZ(num1,K,design_num)+ones(1,design_num);
            for l=1:1:L
                temp1=v(j,l);
                temp2=0;
                for i=1:1:design_num
                    temp1=temp1*z(i,temp_designcluster(i));
                    if (i~=best_design(j))
                        temp3=(mu(best_design(j),j,temp_designcluster(best_design(j)),l)-mu(i,j,temp_designcluster(i),l))^2/(sigma(best_design(j),j,temp_designcluster(best_design(j)),l)^2+sigma(i,j,temp_designcluster(i),l)^2);
                        if temp2==0
                            temp2=temp3;
                        else
                            temp2=min(temp2,temp3);
                        end
                    end
                end
                temp1=temp1*temp2;
                temp=temp+temp1;
            end
        end
        if j==1
            VFA=temp;
        else
            VFA=min(VFA,temp);
        end
    end
end

CS=zeros(1,context_num);
for j=1:1:context_num
    if best_design(j)==true_best(j)
        CS(j)=1;
    end
end