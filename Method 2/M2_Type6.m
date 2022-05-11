
%Prop_Method = 2*Hmax+2+R1+2*R; %First phase 2t slots, 2 for broadcasts by BS, R1 second phase, R third Phase
clc
clear
N_limit=50;
Max_iter = 500;
Prop_Method = zeros(1,Max_iter);
T=6;%Types of users
Type = zeros(T,N_limit);
Nt = zeros(1,T);
Exp_time_slots = zeros(1,10);
b = 6; % slot width
ncc = 6; % Number of case  which require additional slots in next phase
%%%% One of (2,3),  One of (5,6)  ,{1, One of (4,5,6)},  {4, One of (1,2,3)},  1,  4
len = ceil(log2((ncc) + 2));
load 4dev.txt %Load values from 2 types of nodes' result
Four_Dev_Result = X4dev; %Load command saves the values in X4dev by default.
%q = zeros(1,T);
% Active_N = zeros(1,T);
%n = zeros(1,T);
for t=1:T
    Type(t,1:N_limit)=1:N_limit; %Set of T type of nodes
end
for t=1:T
    Nt(t) = length(Type(t,1:N_limit)); %total t^th type of nodes
end

Hmax = ceil(log2(max(Nt))); % Maximum Hash Value.

r=0;
for q=0.02:0.02:0.2
    r=r+1;
    for iter =1:Max_iter
        temp=zeros(Hmax,T);
        for t=1:T
            % temp(1:Hmax,t) = Active_Nodes(Nt(t),q(t),temp(1:Hmax,t)); %Set of hash values for each type of user.
            temp(1:Hmax,t) = Active_Nodes(Nt(t),q,temp(1:Hmax,t)); %Set of hash values for each type of user.
            % n(t) = length(Active_N(t));
        end
        
        R1=0;
        R2 = 0;
        T4C = 0;
        
        for i=1:Hmax
            % 3 Collision Case      
            if ( temp(i,3)>=2 || temp(i,6)>=2 || (temp(i,3)==1 && temp(i,6)==1) )
                    T4C =T4C + 1;
            else %If the above condition is true, 3C always occurs.

                if ( (temp(i,3)==1 && temp(i,6)==0) || (temp(i,3)==0 && temp(i,6)==1) )
                    if( (temp(i,1) + temp(i,2))>=1 && (temp(i,4) + temp(i,5))>=1 && (temp(i,2) + temp(i,5))>=1 )
                        T4C =T4C + 1;
                    end
                    
                elseif (temp(i,3)==0 && temp(i,6)==0)
                    if( (temp(i,1) + temp(i,2))>=2 && (temp(i,4) + temp(i,5))>=2 && (temp(i,2) + temp(i,5)>=2) )
                        T4C =T4C + 1;
                    end
                 end
           % else
            
                switch 1

                     %%%% Two C cases
                    case ((((temp(i,2) >= 2)&&(temp(i,3) == 0))) && (temp(i,4) == 0) && (temp(i,5) == 0) && (temp(i,6) == 0)) % (C, C, 0) case, Type 1 can be anything
                        R1=R1+1;
                    case ((((temp(i,2) >= 1)&&(temp(i,3) == 1))) && (temp(i,4) == 0) && (temp(i,5) == 0) && (temp(i,6) == 0)) % (C, C, alpha) case, Type 1 can be anything
                        R1=R1+1;  
                    case ((((temp(i,2) >= 1)&&(temp(i,3) == 0))) && (temp(i,4) == 0) && (temp(i,5) == 0) && (temp(i,6) == 1)) % (C, C, beta) case, Type 1 can be anything
                        R2=R2+1;
                    case ((((temp(i,2) >= 2)&&(temp(i,3) == 0))) && (temp(i,4) == 1) && (temp(i,5) == 0) && (temp(i,6) == 0)) % (C, C, beta) case, Type 1 can be anything
                        R2=R2+1;
                    case ((((temp(i,2) >= 2)&&(temp(i,3) == 0))) && (temp(i,4) == 0) && (temp(i,5) == 1) && (temp(i,6) == 0)) % (C, C, beta) case, Type 1 can be anything
                        R2=R2+1;
                    case ((temp(i,1) >= 1) && ((temp(i,2) == 1) && (temp(i,3) == 0)) && (temp(i,4) == 0) && (temp(i,5) == 1) && (temp(i,6) == 0)) % (C, C, beta) case
                        R2=R2+1;
                    case ((temp(i,1) >= 1) &&(((temp(i,2) == 1)&&(temp(i,3) == 0))) && (temp(i,4) >= 2) && (temp(i,5) == 0) && (temp(i,6) == 0)) % (C, alpha, C) case.
                        R1=R1+1;
                    case ((temp(i,1) >= 1) &&(((temp(i,2) == 0)&&(temp(i,3) == 1))) && (temp(i,4) >= 1) && (temp(i,5) == 0) && (temp(i,6) == 0)) % (C, alpha, C) case.
                        R1=R1+1;
                    case ((temp(i,1) >= 2) &&(((temp(i,2) == 0)&&(temp(i,3) == 0))) && (temp(i,4) >= 1) && (temp(i,5) == 1) && (temp(i,6) == 0)) % (C, beta, C) case.
                        R1=R1+1;
                    case ((temp(i,1) >= 1) &&(((temp(i,2) == 0)&&(temp(i,3) == 0))) && (temp(i,4) >= 1) && (temp(i,5) == 0) && (temp(i,6) == 1)) % (C, beta, C) case.
                        R1=R1+1;
                    case ((((temp(i,5) >= 2)&&(temp(i,3) == 0))) && (temp(i,1) == 0) && (temp(i,2) == 0) && (temp(i,6) == 0)) % (0, C, C) case, Type 4 can be anything
                        R1=R1+1;
                    case ((((temp(i,5) >= 1)&&(temp(i,3) == 1))) && (temp(i,1) == 0) && (temp(i,2) == 0) && (temp(i,6) == 0)) % (alpha, C, C) case, Type 4 can be anything
                        R2=R2+1;
                    case ((((temp(i,5) >= 2)&&(temp(i,3) == 0))) && (temp(i,1) == 1) && (temp(i,2) == 0) && (temp(i,6) == 0)) % (alpha, C, C) case, Type 4 can be anything
                        R2=R2+1;
                    case ((temp(i,4) >= 1) && (((temp(i,5) == 1)&&(temp(i,3) == 0))) && (temp(i,2) == 1) && (temp(i,1) == 0) && (temp(i,6) == 0)) % (alpha, C, C) case
                        R2=R2+1;
                    case ((((temp(i,5) >= 2)&&(temp(i,3) == 0))) && (temp(i,2) == 1) && (temp(i,1) == 0) && (temp(i,6) == 0)) % (alpha, C, C) case
                        R2=R2+1;
                        
                    case ((((temp(i,5) >= 1)&&(temp(i,6) == 1))) && (temp(i,1) == 0) && (temp(i,2) == 0) && (temp(i,3) == 0)) % (beta, C, C) case, Type 4 can be anything
                        R1=R1+1;

                   %%%% One C cases
                    case ((temp(i,1) == 0) && (((temp(i,2) == 1)&&(temp(i,3) == 0))) && (temp(i,4) >= 2) && (temp(i,5) == 0) && (temp(i,6) == 0)) % (alpha,alpha, C) case
                        R1=R1+1;
                    case ((temp(i,1) == 0) && (((temp(i,2) == 0)&&(temp(i,3) == 1))) && (temp(i,4) >= 1) && (temp(i,5) == 0) && (temp(i,6) == 0)) % (alpha,alpha, C) case
                        R1=R1+1;
                    case ((temp(i,1) >= 2) && (((temp(i,2) == 0)&&(temp(i,3) == 0))) && (temp(i,4) == 0) && (temp(i,5) == 1) && (temp(i,6) == 0)) % (C, beta, beta) case
                        R1=R1+1;
                    case ((temp(i,1) >= 1) && (((temp(i,2) == 0)&&(temp(i,3) == 0))) && (temp(i,4) == 0) && (temp(i,5) == 0) && (temp(i,6) == 1)) % (C, beta, beta) case
                        R1=R1+1;
                end %switch end
            end %if end
        end %for end
        
        Prop_Method(iter) = 3*Hmax+ (len*Hmax)/b + R1 + 2*R2 + + T4C*(2 + (Four_Dev_Result(r)/Hmax)); %2 slots for Type 3 and Type 6 , /Hmax bcoz in prev case, results calculated for Hmax blocks.      
        display(iter)
    end %iter end
    display(q)
    Exp_time_slots(r) = mean(Prop_Method);
end %q end
%Exp_time_slots = mean(Prop_Method);
T_Rep(1:10) = T*Hmax;
%plot(1:10,Exp_time_slots,Three_Rep)
x=.02:.02:.2; %case of N variation
plot(x,Exp_time_slots,'-k',x,T_Rep,'--k')%,x,u,':.k')
set(gca,'FontSize',18,'FontName','Times New Roman')
title(sprintf('Number of slots required for T = %d',T),'FontSize',20)
xlabel('q','FontSize',20)%Probability with which a node is active
ylabel('Number of slots required','FontSize',20)
sum1 = 'Proposed Method';
sum2 = 'T Repetitions of LoF';
Legend = cell(2,1);
Legend{1,1} = (sprintf('%s',sum1));
Legend{2,1} = (sprintf('%s', sum2));
legend(Legend,'FontSize',18)
savefig('Slots_Required_6');
save 6dev.txt Exp_time_slots -ascii %To save the values, which can be used in 8dev, 10dev etc.
%display(Exp_time_slots);
%display(Three_Rep);
