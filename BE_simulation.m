function [next_state]=BE_simulation(current_state,design,design_para,context)

%% context setting
mortality=1/(85-context(1))/12; %context(1) age
doable=1-(context(1)-45)*0.00225;
%% design1
    if design==1
        %% context setting
        effect=0.5+(design_para(1)-75)*0.003-(context(2)-120)*0.005;
        complication=0.025+(design_para(1)-75)*0.0005-(context(2)-120)*0.001;
    
        %% Markov
        Tran_matrix=zeros(8,8);
        % 1:Drug 2:Complication 3:BE 4:Cancer 5:Surgery 6:Inoperable 7:Death 8:Post 
        Tran_matrix(1,7)=mortality;  % all-cause mortality
        Tran_matrix(1,4)=0.005*(1-Tran_matrix(1,7))*(1-effect);
        Tran_matrix(1,2)=complication*(1-Tran_matrix(1,7));
        Tran_matrix(1,1)=1-Tran_matrix(1,2)-Tran_matrix(1,4)-Tran_matrix(1,7);

        Tran_matrix(2,7)=0.025; % complication mortality
        Tran_matrix(2,3)=1-Tran_matrix(2,7); 

        Tran_matrix(3,7)=mortality;  % all-cause mortality
        Tran_matrix(3,4)=0.005*(1-Tran_matrix(3,7));
        Tran_matrix(3,3)=1-Tran_matrix(3,4)-Tran_matrix(3,7);

        Tran_matrix(4,7)=mortality;  % all-cause mortality
        Tran_matrix(4,5)=doable*(1-Tran_matrix(4,7));
        Tran_matrix(4,6)=1-Tran_matrix(4,5)-Tran_matrix(4,7);

        Tran_matrix(5,7)=0.02;   %resection mortality
        Tran_matrix(5,8)=0.8;    %recovery
        Tran_matrix(5,6)=1-Tran_matrix(5,7)-Tran_matrix(5,8);

        Tran_matrix(6,7)=0.1;    %cancer mortality
        Tran_matrix(6,6)=1-Tran_matrix(6,7);

        Tran_matrix(8,7)=mortality;  % all-cause mortality
        Tran_matrix(8,8)=1-Tran_matrix(8,7);

        %%one step
        temp=rand(1);
        if current_state==1
           if temp<=Tran_matrix(1,7)
               next_state=7;
           elseif temp<=Tran_matrix(1,7)+Tran_matrix(1,4)
               next_state=4;
           elseif temp<=Tran_matrix(1,7)+Tran_matrix(1,4)+Tran_matrix(1,2)
               next_state=2;
           else
               next_state=1;
           end
        elseif current_state==2
           if temp<=Tran_matrix(2,7)
               next_state=7;
           else
               next_state=3;
           end        
        elseif current_state==3
           if temp<=Tran_matrix(3,7)
               next_state=7;
           elseif temp<=Tran_matrix(3,7)+Tran_matrix(3,4)
               next_state=4;
           else
               next_state=3;
           end 
        elseif current_state==4
           if temp<=Tran_matrix(4,7)
               next_state=7;
           elseif temp<=Tran_matrix(4,7)+Tran_matrix(4,5)
               next_state=5;
           else
               next_state=6;
           end
        elseif current_state==5
           if temp<=Tran_matrix(5,7)
               next_state=7;
           elseif temp<=Tran_matrix(5,7)+Tran_matrix(5,8)
               next_state=8;
           else
               next_state=6;
           end
        elseif current_state==6
           if temp<=Tran_matrix(6,7)
               next_state=7; 
           else
               next_state=6;
           end
        elseif current_state==8
           if temp<=Tran_matrix(8,7)
               next_state=7; 
           else
               next_state=8;
           end
        else
           next_state=7; 
        end
    end
%% design2
    if design==2
        %% context setting
        effect=0.5+(design_para(2)-9)*0.0417-(context(2)-120)*0.0025;
        complication=0.04+(design_para(2)-9)*0.01-(context(2)-120)*0.001;
        
        %% Markov
        Tran_matrix=zeros(8,8);
        % 1:Drug 2:Complication 3:BE 4:Cancer 5:Surgery 6:Inoperable 7:Death 8:Post
        Tran_matrix(1,7)=mortality;  % all-cause mortality
        Tran_matrix(1,4)=0.005*(1-Tran_matrix(1,7))*(1-effect);
        Tran_matrix(1,2)=complication*(1-Tran_matrix(1,7));
        Tran_matrix(1,1)=1-Tran_matrix(1,2)-Tran_matrix(1,4)-Tran_matrix(1,7);
        
        Tran_matrix(2,7)=0.025; % complication mortality
        Tran_matrix(2,3)=1-Tran_matrix(2,7);
        
        Tran_matrix(3,7)=mortality;  % all-cause mortality
        Tran_matrix(3,4)=0.005*(1-Tran_matrix(3,7));
        Tran_matrix(3,3)=1-Tran_matrix(3,4)-Tran_matrix(3,7);
        
        Tran_matrix(4,7)=mortality;  % all-cause mortality
        Tran_matrix(4,5)=doable*(1-Tran_matrix(4,7));
        Tran_matrix(4,6)=1-Tran_matrix(4,5)-Tran_matrix(4,7);
        
        Tran_matrix(5,7)=0.02;   %resection mortality
        Tran_matrix(5,8)=0.8;    %recovery
        Tran_matrix(5,6)=1-Tran_matrix(5,7)-Tran_matrix(5,8);
        
        Tran_matrix(6,7)=0.1;    %cancer mortality
        Tran_matrix(6,6)=1-Tran_matrix(6,7);
        
        Tran_matrix(8,7)=mortality;  % all-cause mortality
        Tran_matrix(8,8)=1-Tran_matrix(8,7);
        
        %%one step
        temp=rand(1);
        if current_state==1
            if temp<=Tran_matrix(1,7)
                next_state=7;
            elseif temp<=Tran_matrix(1,7)+Tran_matrix(1,4)
                next_state=4;
            elseif temp<=Tran_matrix(1,7)+Tran_matrix(1,4)+Tran_matrix(1,2)
                next_state=2;
            else
                next_state=1;
            end
        elseif current_state==2
            if temp<=Tran_matrix(2,7)
                next_state=7;
            else
                next_state=3;
            end
        elseif current_state==3
            if temp<=Tran_matrix(3,7)
                next_state=7;
            elseif temp<=Tran_matrix(3,7)+Tran_matrix(3,4)
                next_state=4;
            else
                next_state=3;
            end
        elseif current_state==4
            if temp<=Tran_matrix(4,7)
                next_state=7;
            elseif temp<=Tran_matrix(4,7)+Tran_matrix(4,5)
                next_state=5;
            else
                next_state=6;
            end
        elseif current_state==5
            if temp<=Tran_matrix(5,7)
                next_state=7;
            elseif temp<=Tran_matrix(5,7)+Tran_matrix(5,8)
                next_state=8;
            else
                next_state=6;
            end
        elseif current_state==6
            if temp<=Tran_matrix(6,7)
                next_state=7;
            else
                next_state=6;
            end
        elseif current_state==8
            if temp<=Tran_matrix(8,7)
                next_state=7;
            else
                next_state=8;
            end
        else
            next_state=7;
        end
    end
end
