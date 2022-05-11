function [a, b, c, d, e, f, g, h]= bw2(e_hat, p_hat, n_hat, Ea, Pa, Na, Sset, M, Stack,Trem)
%Simulates the Broadcast window 2
Se=[];
Sp=[];
Sn=[];
We=3;
Wp=2;
Wn=1;
%Pe=0;
%Pp=0;
%Pn=0;
CM1=cell(M,1);
N=length(Sset);

temp=round((e_hat*We*N)/(e_hat*We + p_hat*Wp + n_hat*Wn));         %Function, to be decided later
for i=1:temp
    Se=union(Se,Sset(i));                         %Channels for emergency nodes
end
Sset=setdiff(Sset,Se);

temp=round((p_hat*Wp*N)/(e_hat*We + p_hat*Wp + n_hat*Wn));         %Function, to be decided later
for i=1:temp
    Sp=union(Sp,Sset(i));                          %Channels for periodic nodes
end
Sn=setdiff(Sset,Sp);                            %Channels for normal nodes

l1=length(Se);                  
l2=length(Sp);                   
l3=length(Sn); 

Pe=l1/e_hat;                   %emergency nodes attempt probability
Pe=min(Pe,0.8);
Pp=l2/p_hat;                   %periodic nodes attempt probability 
Pp=min(Pp,0.8);
Pn=l3/n_hat;                   %normal nodes attempt probability
Pn=min(Pn,0.8);

for i=1:length(Ea)   
    if(l1==0)
    else
    ran=randi(l1,1);            %Choose a channel from Se at random
    temp_ch=Se(ran);            %Allocate the channel to Ea(i)
    CM1{temp_ch}=[CM1{temp_ch} Ea(i)];  %Update CM
    end
end
    
for i=1:length(Pa)
    if(l2==0)
    else
    ran=randi(l2,1);           %Choose a channel from Sp at random
    temp_ch=Sp(ran);           %Allocate the channel to Pa(i)
    CM1{temp_ch}=[CM1{temp_ch} Pa(i)]; %Update CM
    end
end
    
for i=1:length(Na)
    if(l3==0)
    else
    ran=randi(l3,1);           %Choose a channel from Sn at random
    temp_ch=Sn(ran);           %%Allocate the channel to Na(i)
    CM1{temp_ch}=[CM1{temp_ch} Na(i)]; %Update CM
    end
end

%to move periodic nodes from a busy channel-- check the code
Nonempty=find(~cellfun(@isempty,Stack));
for i=1:length(Nonempty)
     if(~ismember(Nonempty(i),Sp))
         if(l2==0)
            Stack(Nonempty(i))={[]};
         else
             ran=randi(l2,1);           %Choose a channel from Sp at random
             temp_ch=Sp(ran);
             if(length(Stack{temp_ch})+ length(Stack{Nonempty(i)})<2*Trem)
                 Stack{temp_ch}=[Stack{temp_ch} Stack{Nonempty(i)}];
                 %%Stack(temp_ch,2)=Stack(Nonempty(i),2);   Stack changed to 1D cell later
             end
             Stack(Nonempty(i))={[]};
         end
     end
end
    

a=Se;
b=Sp;
c=Sn;
d=Pe;
e=Pp;
f=Pn;
g=CM1;
h=Stack;