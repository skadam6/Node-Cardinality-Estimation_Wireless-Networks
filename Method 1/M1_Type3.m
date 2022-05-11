
%Prop_Method = 2*Hmax+2+K+2*R; %First phase 2t slots, 2 for broadcasts by BS, K second phase, R third Phase
clc
clear 
N_limit=50;
Max_iter = 500;
Prop_Method = zeros(1,Max_iter);
T=3;%Types of users
Type = zeros(T,N_limit);
Nt = zeros(1,T);
Exp_time_slots = zeros(1,10);
q_max = .1;
%q = zeros(1,T);
% Active_N = zeros(1,T);
n = zeros(1,T);
for t=1:T
    Type(t,1:N_limit)=1:N_limit; %Set of T type of nodes
end
for t=1:T
    Nt(t) = length(Type(t,1:N_limit)); %total t^th type of nodes
end

Hmax = ceil(log2(max(Nt))); % Maximum Hash Value.
% for t=1:T
%     q(t) = 0.05;%its a case. In future different prob values can be used.
% end
r=0;
for q=0.01:0.01:q_max
    r=r+1;
for iter =1:Max_iter
temp=zeros(Hmax,T);
for t=1:T    
   % temp(1:Hmax,t) = Active_Nodes(Nt(t),q(t),temp(1:Hmax,t)); %Set of hash values for each type of user.
    temp(1:Hmax,t) = Active_Nodes(Nt(t),q,temp(1:Hmax,t)); %Set of hash values for each type of user.
   % n(t) = length(Active_N(t));
end

K=0;
R=0;
for i=1:Hmax
        %Considered all cases where collision might happen
        if (((temp(i,1)==1)&&(sum(temp(i,2:T)>=1) == (T-1)))||((temp(i,1)==0)&&(sum(temp(i,2:T)>=2) == (T-1)))||((temp(i,1)>=2)))            
            K=K+1;          
        end
        if ((temp(i,1)>=2))
            R=R+1;
        end
end

%Prop_Method(iter) = (T-1)*Hmax+2+K+(T-1)*R; %First phase 2t slots, 2 for broadcasts by BS, K second phase, R third Phase
if (K == 0)
    Prop_Method(iter) = (T-1)*Hmax+1; %First phase 2t slots, 1 for BP after 1st phase by BS, K second phase, 1 for BP after 2nd phase(it may be zero). R third Phase
else
    Prop_Method(iter) = (T-1)*Hmax+1+K+1+(T-1)*R; %First phase 2t slots, 1 for BP after 1st phase by BS, K second phase, 1 for BP after 2nd phase(it may be zero). R third Phase
end
display(iter)
end
display(q)
Exp_time_slots(r) = mean(Prop_Method);
end
%Exp_time_slots = mean(Prop_Method);
T_Rep(1:10) = T*Hmax;
%plot(1:10,Exp_time_slots,Three_Rep)
x=.01:.01:q_max; %case of N variation
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
savefig('Slots_Required_3_Method2');
save 3dev.txt Exp_time_slots -ascii
%display(Exp_time_slots);
%display(Three_Rep);
