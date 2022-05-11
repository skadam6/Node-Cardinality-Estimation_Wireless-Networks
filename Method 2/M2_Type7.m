
%Prop_Method = 2*Hmax+2+K+2*R; %First phase 2t slots, 2 for broadcasts by BS, K second phase, R third Phase
clc
clear
N_limit=50;
Max_iter = 500;
Prop_Method = zeros(1,Max_iter);
T=7;%Types of users
Type = zeros(T,N_limit);
Nt = zeros(1,T);
Exp_time_slots = zeros(1,10);
b = 6; % slot width
ncc = 11; % Number of case  which require additional slots in next phase
len = ceil(log2((ncc) + 2));
load 6dev.txt %Load values from 2 types of nodes' result
Six_Dev_Result = X6dev; %Load command saves the values in X6dev by default.
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
            temp(1:Hmax,t) = Active_Nodes(Nt(t),q,temp(1:Hmax,t)); %Set of hash values for each type of user.
        end
        
        K=0;
        R=0;
        R1 = 0;
        R2 = 0;
        R3 = 0;
		R4 = 0;
        T4C = 0;
        
        for i=1:Hmax
            %%%% 4 Collision Case
            % ======================================================= 4 Collision=============================================================
          
            if ( temp(i,4)>=2 )
                T4C = T4C +1;
            
            else
                if ( (temp(i,4)==1) )
                    if( sum(temp(i,1:3))>=1 && sum(temp(i,5:7))>=1 && sum(temp(i,[2,3,7]))>=1 && sum(temp(i,[3,6,7]))>=1 )
                        T4C = T4C +1;
                    end

                elseif ( (temp(i,4)==0) )
                    if( sum(temp(i,1:3))>=2 && sum(temp(i,5:7))>=2 && sum(temp(i,[2,3,7]))>=2 && sum(temp(i,[3,6,7]))>=2 )
                       T4C = T4C +1;
                    end
                end

                switch 1

                    %%%% 3 C cases
                    % =======================================================3 Collision=============================================================
                    case ( (temp(i,3) >= 2) && (temp(i,4)==0) && (temp(i,5)==0) && (temp(i,6)==0) && (temp(i,7)==0)) % CCC0 Case, 2 slots for Type 1 and 2
                        R2 = R2+1;

                    case ( (temp(i,3) >= 1) && (temp(i,4)==1) && (temp(i,5)==0) && (temp(i,6)==0) && (temp(i,7)==0)) % CCCa Case, 2 slots for Type 1 and 2                  
                        R2 = R2+1;

                    case ((temp(i,3) >= 2) && (temp(i,4)==0) && (temp(i,5)==1) && (temp(i,6)==0) && (temp(i,7)==0))   % CCCb case, 3 slots for Type (1; 2; One among {5; 6; 7}), 5==1
                        R3 = R3+1;
                    case ((temp(i,3) >= 2) && (temp(i,4)==0) && (temp(i,5)==0) && (temp(i,6)==1) && (temp(i,7)==0))   % CCCb 6 == 1
                        R3 = R3+1;
                    case ((temp(i,2) >= 1) && (temp(i,3) == 1) && (temp(i,4)==0) && (temp(i,5)==0) && (temp(i,6)==1) && (temp(i,7)==0)) % CCCb 6 == 1
                        R3 = R3+1;
                    case ((temp(i,3) >= 2) && (temp(i,4)==0) && (temp(i,5)==0) && (temp(i,6)==0) && (temp(i,7)==1)) % CCCb  7 == 1
                        R3 = R3+1;    
                    case ( ((temp(i,1) >= 1) || (temp(i,2) >= 1)) && (temp(i,3) == 1) && (temp(i,4)==0) && (temp(i,5)==0) && (temp(i,6)==0) && (temp(i,7)==1)) % CCCb 7 == 1
                        R3 = R3+1;    
                    

                    case ((temp(i,2) >= 2) && (temp(i,3)==0) && (temp(i,4)==0) && (temp(i,5)>= 2) && (temp(i,6)==0) && (temp(i,7)==0))   % CC0C, 1 slot for Type 1 
                        R1 = R1+1;

                    case ((temp(i,2) >= 1) && (temp(i,3)==1) && (temp(i,4)==0) && (temp(i,5)>= 2) && (temp(i,6)==0) && (temp(i,7)==0))   % CCaC, 2 slot for Type 1 and one among 3,4 
                        R2 = R2+1;
                    case ((temp(i,2) >= 1) && (temp(i,3)==0) && (temp(i,4)==1) && (temp(i,5)>= 1) && (temp(i,6)==0) && (temp(i,7)==0))   % CCaC, 2 slot for Type 1 and one among 3,4 
                        R2 = R2+1;   

                    case ((temp(i,2) >= 2) && (temp(i,3)==0) && (temp(i,4)==0) && (temp(i,5)>= 1) && (temp(i,6)==1) && (temp(i,7)==0))   % CCbC, 2 slot for Type 1 and One among 6;7
                        R2 = R2+1;   
                    case ((temp(i,2) >= 2) && (temp(i,3)==0) && (temp(i,4)==0) && (temp(i,5)>= 1) && (temp(i,6)==0) && (temp(i,7)==1))   % CCbC, 2 slot for Type 1 and One among 6;7
                        R2 = R2+1;    
                    case ((temp(i,1) >= 1) && (temp(i,2) == 1) && (temp(i,3)==0) && (temp(i,4)==0) && (temp(i,5)>= 1) && (temp(i,6)==0) && (temp(i,7)==1))   % CCbC, 2 slot for Type 1 and One among 6;7
                        R2 = R2+1;     
                         

                    case ((temp(i,6) >= 2) && (temp(i,7)==0) && (temp(i,1)>= 2) && (temp(i,2)==0) && (temp(i,3)==0) && (temp(i,4)==0))   % C0CC, 1 slot for Type 5 
                        R1 = R1+1;                           
                                                             
                    case ((temp(i,6) >= 2) && (temp(i,7)==0) && (temp(i,1)>= 1) && (temp(i,2)==1) && (temp(i,3)==0) && (temp(i,4)==0))   % CaCC, 2 slot for Type 5 and One among 2,3,4
                        R2 = R2+1;                           
                    case ((temp(i,6) >= 2) && (temp(i,7)==0) && (temp(i,1)>= 1) && (temp(i,2)==0) && (temp(i,3)==1) && (temp(i,4)==0))   % CaCC, 2 slot for Type 5 and One among 2,3,4
                        R2 = R2+1;    
                    case ((temp(i,5) >= 1) && (temp(i,6) == 1) && (temp(i,7)==0) && (temp(i,1)>= 1) && (temp(i,2)==0) && (temp(i,3)==1) && (temp(i,4)==0))   % CaCC, 2 slot for Type 5 and One among 2,3,4
                        R2 = R2+1;     
                    case ((temp(i,6) >= 1) && (temp(i,7)==0) && (temp(i,1)>= 1) && (temp(i,2)==0) && (temp(i,3)==0) && (temp(i,4)==1))   % CaCC, 2 slot for Type 5 and One among 2,3,4
                        R2 = R2+1;     

                    case ((temp(i,6) >= 1) && (temp(i,7)==1) && (temp(i,1)>= 2) && (temp(i,2)==0) && (temp(i,3)==0) && (temp(i,4)==0))   % CbCC, 1 slot for Type 5
                        R1 = R1+1;
                        

                    case ( (temp(i,7) >= 2)&& (temp(i,1)==0) && (temp(i,2)==0) && (temp(i,3)==0) && (temp(i,4)==0) ) % 0CCC Case, 2 slots for Type 5 and 6
                        R2 = R2+1;   

                    case ((temp(i,7) >= 2) && (temp(i,1)==1) && (temp(i,2)==0) && (temp(i,3)==0) && (temp(i,4)==0))   % aCCC case, 4 slots for Type (5; 6; One among {1;2;3;4}), 1==1
                        R4 = R4+1;         
                    case ((temp(i,7) >= 2) && (temp(i,1)==0) && (temp(i,2)==1) && (temp(i,3)==0) && (temp(i,4)==0))   % aCCC case
                        R4 = R4+1;
                    case ((temp(i,6) >= 1) && (temp(i,7) == 1) && (temp(i,1)==0) && (temp(i,2)==1) && (temp(i,3)==0) && (temp(i,4)==0)) % aCCC case,
                        R4 = R4+1;
                    case ((temp(i,7) >= 2) && (temp(i,1)==0) && (temp(i,2)==0) && (temp(i,3)==1) && (temp(i,4)==0)) % aCCC case,
                        R4 = R4+1;    
                    case ( ((temp(i,5) >= 1) || (temp(i,6) >= 1)) && (temp(i,7) == 1) && (temp(i,1)==0) && (temp(i,2)==0) && (temp(i,3)==1) && (temp(i,4)==0)) % aCCC case,
                        R4 = R4+1;    
                    case ((temp(i,7) >= 1) && (temp(i,1)==0) && (temp(i,2)==0) && (temp(i,3)==0) && (temp(i,4)==1)) % aCCC case,
                        R4 = R4+1;   

                    
                   %%%% Two C cases
                   % =====================================================2 Collision=============================================================
                   case ((temp(i,2) >= 2) && (temp(i,3)==0) && (temp(i,4)==0) && (temp(i,5)== 0) && (temp(i,6)==0) && (temp(i,7)==0))   % CC00, 1 slot for Type 1
                        R1 = R1+1;
                   case ((temp(i,2) >= 2) && (temp(i,3)==0) && (temp(i,4)==0) && (temp(i,5)== 1) && (temp(i,6)==0) && (temp(i,7)==0))   % CC0b, 1 slot for Type 1
                        R1 = R1+1; 
                   case ((temp(i,2) >= 1) && (temp(i,3)==1) && (temp(i,4)==0) && (temp(i,5)== 0) && (temp(i,6)==0) && (temp(i,7)==0))   % CCa0, 1 slot for Type 1
                        R1 = R1+1; 
                   case ((temp(i,2) >= 1) && (temp(i,3)==0) && (temp(i,4)==1) && (temp(i,5)== 0) && (temp(i,6)==0) && (temp(i,7)==0))   % CCaa, 1 slot for Type 1
                        R1 = R1+1;      
                   case ((temp(i,2) >= 1) && (temp(i,3)==1) && (temp(i,4)==0) && (temp(i,5)== 1) && (temp(i,6)==0) && (temp(i,7)==0))   % CCab, 1 slot for Type 1
                        R1 = R1+1;     

                   case ((temp(i,2) >= 2) && (temp(i,3)==0) && (temp(i,4)==0) && (temp(i,5)== 0) && (temp(i,6)==1) && (temp(i,7)==0) )   % CCbb, 2 slot for Type 1 ,One among 6,7
                        R2 = R2+1;
                   case ((temp(i,2) >= 2) && (temp(i,3)==0) && (temp(i,4)==0) && (temp(i,5)== 0) && (temp(i,6)==0) && (temp(i,7)==1) )   % CCbb, 2 slot for Type 1 ,One among 6,7
                        R2 = R2+1;    
                   case ((temp(i,1) >= 1) && (temp(i,2) == 1) && (temp(i,3)==0) && (temp(i,4)==0) && (temp(i,5)== 0) && (temp(i,6)==0) && (temp(i,7)==1))   % CCbb, 2 slot for Type 1 ,One among 6,7
                        R2 = R2+1;     
                           

                   case ((temp(i,1) >= 1) && (temp(i,2) == 0) && (temp(i,3)==1) && (temp(i,4)==0) && (temp(i,5)>= 2) && (temp(i,6)==0) && (temp(i,7)==0))   % CaaC, 1 slot for One among 3,4
                        R1 = R1+1;
                   case ((temp(i,1) >= 1) && (temp(i,2) == 0) && (temp(i,3)==0) && (temp(i,4)==1) && (temp(i,5)>= 1) && (temp(i,6)==0) && (temp(i,7)==0))   % CaaC, 1 slot for One among 3,4
                        R1 = R1+1;     

                        

                   case ((temp(i,6) >= 2) && (temp(i,7)==0) && (temp(i,1)== 0) && (temp(i,2)==0) && (temp(i,3)==0) && (temp(i,4)==0))   % 00CC, 1 slot for Type 5
                        R1 = R1+1;
                   case ((temp(i,6) >= 1) && (temp(i,7)==1) && (temp(i,1)== 0) && (temp(i,2)==0) && (temp(i,3)==0) && (temp(i,4)==0))   % 0bCC, 1 slot for Type 5
                        R1 = R1+1;     
                   case ((temp(i,6) >= 2) && (temp(i,7)==0) && (temp(i,1)== 1) && (temp(i,2)==0) && (temp(i,3)==0) && (temp(i,4)==0))   % a0CC, 1 slot for Type 5
                        R1 = R1+1;    

                   case ((temp(i,6) >= 2) && (temp(i,7)==0) && (temp(i,1)== 0) && (temp(i,2)==1) && (temp(i,3)==0) && (temp(i,4)==0))   % aaCC, 2 slot for Type 1 ,One among 2,3,4
                        R2 = R2+1;
                   case ((temp(i,6) >= 2) && (temp(i,7)==0) && (temp(i,1)== 0) && (temp(i,2)==0) && (temp(i,3)==1) && (temp(i,4)==0))   % aaCC, 2 slot for Type 1 ,One among 2,3,4
                        R2 = R2+1;    
                   case ((temp(i,5) >= 1) && (temp(i,6) == 1) && (temp(i,7)==0) && (temp(i,1)== 0) && (temp(i,2)==0) && (temp(i,3)==1) && (temp(i,4)==0))  % aaCC, 2 slot for Type 1 ,One among 2,3,4
                        R2 = R2+1;     
                   case ((temp(i,6) >= 1) && (temp(i,7)==0) && (temp(i,1)== 0) && (temp(i,2)==0) && (temp(i,3)==0) && (temp(i,4)==1))   % aaCC, 2 slot for Type 1 ,One among 2,3,4
                        R2 = R2+1;     

                   case ((temp(i,6) >= 1) && (temp(i,7)==1) && (temp(i,1)== 1) && (temp(i,2)==0) && (temp(i,3)==0) && (temp(i,4)==0))   % abCC, 1 slot for Type 5
                        R1 = R1+1;   
                      

                   %%%% One C cases
                   % =====================================================2 Collision=============================================================

                    
                    case ((temp(i,5) >= 2) && (temp(i,6) == 0) && (temp(i,7)==0) && (temp(i,1)== 0) && (temp(i,2)==0) && (temp(i,3)==1) && (temp(i,4)==0))  % aaaC, 1 slot for One among 3,4
                        R1 = R1+1;                                               
                    case ((temp(i,5) >= 1) && (temp(i,6) == 0) && (temp(i,7)==0) && (temp(i,1)== 0) && (temp(i,2)==0) && (temp(i,3)==0) && (temp(i,4)==1))  % aaaC, 1 slot for One among 3,4
                        R1 = R1+1;          

                end
            end
        end
   
       Prop_Method(iter) = 4*Hmax+ (len*Hmax)/b + R1 + 2*R2 + 3*R3 + 4*R4 + T4C*(1 + (Six_Dev_Result(r)/Hmax)); %2 slots for Type 4 and Type 8 , /Hmax bcoz in prev case, results calculated for Hmax blocks.      
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
savefig('Slots_Required_7');
save 7dev.txt Exp_time_slots -ascii %To save the values, which can be used in 8dev, 10dev etc.
%display(Exp_time_slots);
%display(Three_Rep);
