
%Prop_Method = 2*Hmax+2+K+2*R; %First phase 2t slots, 2 for broadcasts by BS, K second phase, R third Phase
clc
clear
N_limit=50;
Max_iter = 500;
Prop_Method = zeros(1,Max_iter);
T=5;%Types of users
Type = zeros(T,N_limit);
Nt = zeros(1,T);
Exp_time_slots = zeros(1,10);
b = 6; % slot width
ncc = 6; % Number of case  which require additional slots in next phase
%%%% One of (3,4), One for (1,2)
len = ceil(log2((ncc) + 2));

load 3dev.txt %Load values from 2 types of nodes' result
Three_Dev_Result = X3dev; %Load command saves the values in X2dev by default.
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
        
        
        R1 = 0;
        R2 = 0;
        R3 = 0;
        Rc = 0;
        for i=1:Hmax
            %Considered all cases where collision might happen
            %IN SUM FINCTION ">=" , IT RETURNS THE SUM OF NUMBER OF COLUMNS
            %WHICH SATISFY THAT EQUATION
            switch 1
               % (0, C) case 2 slot required to distnguish between {3,5},
               % first check for 3 then check for 5
               case ((temp(i,1) == 0) && (temp(i,2) == 0) && ((temp(i,3) >= 2)||(temp(i,5) >= 2))  && (temp(i,4) == 0))
                    R2=R2+1;
               case ((temp(i,1) == 0) && (temp(i,2) == 0) && (((temp(i,3) == 1)&&(temp(i,5) > 1))||((temp(i,3) > 1)&&(temp(i,5) == 1))) && (temp(i,4) == 0)) % if 3 = 0 or 1 no extra slot needed
                    R1=R1+1;
               case ((temp(i,1) == 0) && (temp(i,2) == 0) && (temp(i,3) == 1) && (temp(i,4) == 0) && (temp(i,5) == 1))
                    R1=R1+1;
               % (alpha, C) case, in first phase cheach if 1 present or
               % not, if not then two must be there, in second phase cheack
               % for 3 and 5 seprately, therefor 3 slots required in total to resolve this ambiguity 
               case ((temp(i,1) == 1) && (temp(i,2) == 0) && ((temp(i,3) >= 2)||(temp(i,5) >= 2))  && (temp(i,4) == 0))
                    R3=R3+1;
               case ((temp(i,1) == 1) && (temp(i,2) == 0) && (((temp(i,3) == 1)&&(temp(i,5) > 1))||((temp(i,3) > 1)&&(temp(i,5) >= 1))) && (temp(i,4) == 0))
                    R2=R2+1;
               case ((temp(i,1) == 1) && (temp(i,2) == 0) && (temp(i,3) == 1) && (temp(i,4) == 0) && (temp(i,5) == 1))
                    R2=R2+1;
               case ((temp(i,1) == 0) && (temp(i,2) == 1) && ((temp(i,3) >= 1)||(temp(i,5) >= 1)) && (temp(i,4) == 0))
                    R3 = R3 + 1;
               % (beta, C) case   2 slots required to distnguish bwtween {3,5}   
               case ((temp(i,1) == 0) && (temp(i,2) == 0) && ((temp(i,3) >= 1)||(temp(i,5) >= 1)) && (temp(i,4) == 1))
                    R2 = R2 + 1;
               % (C, alpha) case One slot required to distnguish between one among {2,5}    
               case ((temp(i,1) >= 1) && (temp(i,2) == 1) && (temp(i,3) == 0) && (temp(i,4) == 0) && (temp(i,5) == 0))
                    R1=R1+1;
               case ((temp(i,1) >= 2) && (temp(i,2) == 0) && (temp(i,3) == 0) && (temp(i,4) == 0) && (temp(i,5) == 1))
                    R1=R1+1;
               % (C, beta) case One slot required to distnguish between one among {3,4} 
               case ((temp(i,1) >= 1) && (temp(i,2) == 0) && (temp(i,3) == 0) && (temp(i,4) == 1) && (temp(i,5) == 0))
                    R1=R1+1;
               case ((temp(i,1) >= 2) && (temp(i,2) == 0) && (temp(i,3) == 1) && (temp(i,4) == 0) && (temp(i,5) == 0))
                    R1=R1+1;
               % (alpha, alpha) 1 slot for distinguish between (1,5) or 2,
               % if 2 present then (1,5) absent, therfor only 1 slot
               % required
               case ((temp(i,1) == 1) && (temp(i,2) == 0) && (temp(i,3) == 0) && (temp(i,4) == 0) && (temp(i,5) == 1))
                    R1=R1+1;
               case ((temp(i,1) == 1) && (temp(i,2) == 1) && (temp(i,3) == 0) && (temp(i,4) == 0) && (temp(i,5) == 0))
                    R1=R1+1;
               % (C, C) case
               case ((temp(i,2)>=2)) % (C,C) case
                    Rc=Rc+1;
               case ((temp(i,4)>=2)) % (C,C) case
                    Rc=Rc+1;
               case ((temp(i,2) == 1) && (temp(i,4)==1)) % (C,C) case
                    Rc=Rc+1;
               case ( (temp(i,1) >= 1) && (temp(i,2)==1) && ((temp(i,3)>=1)||(temp(i,5)>=1)) && (temp(i,4)==0)) % (C,C) case
                    Rc=Rc+1;
               case ( (temp(i,1) >= 1) && (temp(i,2)==0) && ((temp(i,3)>=1)||(temp(i,5)>=1)) && (temp(i,4)==1)) % (C,C) case
                    Rc=Rc+1;
               case ( (temp(i,1) >= 2) && (temp(i,2)==0) && ((temp(i,3)>=2)||(temp(i,5)>=2)||((temp(i,3)==1)&&(temp(i,5)==1))) && (temp(i,4)==0)) % (C,C) case
                    Rc=Rc+1;
                    
            end
           
        end
        
       Prop_Method(iter) = 2*Hmax+ (len*Hmax)/b + R1+ 2*R2 + 3*R3 +  Rc*(2 + (Three_Dev_Result(r)/Hmax));
        display(iter)
    end
    display(q)
    Exp_time_slots(r) = mean(Prop_Method);
end
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
savefig('Slots_Required_5new');
save 5devnew.txt Exp_time_slots -ascii %To save the values, which can be used in 6dev, 8dev etc.
%display(Exp_time_slots);
%display(Three_Rep);
