clc;
%-------------------- Initialization ---------------------------------% 
%C=[0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1]; %total channels
C=0.15.*ones(1,30);
M=length(C); %total no of Channels
N_limit=50;
alpha=0.2;
E=linspace(1,N_limit,N_limit); %Set of Emergency nodes
P=linspace(1,N_limit,N_limit); %Set of Periodic nodes
Np=length(P); %total periodic nodes
N=linspace(1,N_limit,N_limit); %Set of Normal nodes
Ne=length(E);
Nn=length(N);

%Ne=ceil((1-alpha)*Np); %total emergency nodes
%E=linspace(1,Ne,Ne); %Set of Emergency nodes
%Nn=ceil((1+alpha)*Np); %total normal nodes
%N=linspace(1,Nn,Nn); %Set of Emergency nodes

total_nodes=Ne+Np+Nn;  %total nodes
Pset = []; %Channels used by PUs
Sset = []; %SU usable channel list
Mf=0;  %Size of Sset
mf=0;   %Broadcasting Channel. i.e. first free channel 
lambda_p=1; %Packet arrival rate of periodic nodes
lambda_e=lambda_p*(1-alpha); %Packet arrival rate of emergency nodes
lambda_n=lambda_p*(1+alpha); %Packet arrival rate of normal nodes
Ea=[]; %Set of emergency nodes with active queue
Na=[]; %Set of Normal nodes with active queue
Pa=[]; %Set of periodic nodes with active queue
Hmax= ceil(log2(max(max(Ne,Nn),Np)))+1;% Maximum Hash Value. 
B=zeros(Hmax,3); %Set used in EW
Tframe=50;           %Total frame time
Tsense=0.2;           %Time spent for sensing the channel
Tspent=0;             %Remaining Time
%-------------------- Queues ---------------------------------% 
Qe=zeros(Ne,1);  % Queue for emergency nodes
Qn=zeros(Nn,1);  % Queue for normal nodes
Qp=zeros(Np,2);  % Queue for periodic nodes- structure description in the notebook
Stack=cell(M,1);  %Stack of every channel, needs to be 2D format=(node,frames reserved) 
Prev_set=[];        %Set of periodic nodes who have previosly reserved access 
%-------------------Idealised Protocol Variables-----------------------
Ea_i=[];
Pa_i=[];
Na_i=[];
Qe_i=zeros(Ne,1);
Qp_i=zeros(Np,2);
Qn_i=zeros(Nn,1);
Stack_i=cell(M,1);
Prev_set_i=[];
%----------------------------------------------
run_limit=5;
lambda_set=(0:0.5:5);
size_lambda=length(lambda_set);
F_limit=50;
TP_e=zeros(1,F_limit);
TP_n=zeros(1,F_limit);
TP_p=zeros(1,F_limit);
TP_e_i=zeros(1,F_limit);
TP_n_i=zeros(1,F_limit);
TP_p_i=zeros(1,F_limit);
avgte=zeros(run_limit,size_lambda);
avgtp=zeros(run_limit,size_lambda);
avgtn=zeros(run_limit,size_lambda);
ideal_e=zeros(run_limit,size_lambda);
ideal_p=zeros(run_limit,size_lambda);
ideal_n=zeros(run_limit,size_lambda);

h = waitbar(0,'Initializing..');
h1=waitbar(0,'Arrival rates');
for RUN=1:run_limit
for n=0:0.5:5
    
    lambda_p=n; %Packet arrival rate of periodic nodes
    %lambda_e=lambda_p*(1-alpha); %Packet arrival rate of emergency nodes
    %lambda_n=lambda_p*(1+alpha); %Packet arrival rate of normal nodes
    lambda_e=lambda_p; 
    lambda_n=lambda_p;
    Ea=[]; %Set of emergency nodes with active queue
    Na=[]; %Set of Normal nodes with active queue
    Pa=[]; %Set of periodic nodes with active queue
    Qe=zeros(Ne,1);  % Queue for emergency nodes
    Qn=zeros(Nn,1);  % Queue for normal nodes
    Qp=zeros(Np,2);  % Queue for periodic nodes- structure description in the notebook
    Stack=cell(M,1);  %Stack of every channel, needs to be 2D format=(node,frames reserved) 
    Prev_set=[];        %Set of periodic nodes who have previosly reserved access 
 for frames=1:F_limit
%-------------------- Packet Generation ---------------------------------% 
    for i=1:Ne    %Packet generation for emergency nodes
        k = poissrnd(lambda_e);
        if k ~= 0
            Ea= union(Ea,i);
            Qe(i)= Qe(i)+k;
        end
    end
    for i=1:Nn    %Packet generation for normal nodes
        k = poissrnd(lambda_n);
        if k ~= 0
            Na= union(Na,i);
             Qn(i)= Qn(i)+k;
        end
    end
   for i=1:Np     %Packet generation for periodic nodes
      k = poissrnd(lambda_p);
       if k ~= 0
          if (any(Prev_set==i))
              Qp(i,2)= Qp(i,2)+ k;
            else
                Pa= union(Pa,i);
                Qp(i,1)=Qp(i,1)+k;
          end
       end
    end
      Ea_i=Ea;
      Pa_i=Pa;
      Na_i=Na;
      Qe_i=Qe;
      Qp_i=Qp;
      Qn_i=Qn;
      Prev_set_i=Prev_set;
%-------------------- BW1 ---------------------------------%
    [Pset, Sset]= bw1(C);
    Sset_i=Sset;
    mf=Sset(1);
    Mf=length(Sset);
    Tbw1=Tsense*(mf + 1);  % Time spent in BW1 (+1 is the time taken by BS to sense channel)
    Tbw1=ceil(Tbw1);

%-------------------- EW ---------------------------------%
     [B, e_hat, p_hat, n_hat, E2, E3, phase2, phase3]=EW(Ea,Na,Pa,Hmax,Ne,Np,Nn);
     Tewp1=ceil((ceil(log2(max(max(Ne,Nn),Np)))+1)/Mf);
     x=size(E2);
     Tewp2=ceil(x(1)/Mf);
     x=size(E3);
     Tewp3=ceil(x(1)/Mf);
     Tew=Tewp1+Tewp2+Tewp3;
     Tspent=Tbw1+Tew;
     Trem=Tframe-Tspent;
     Trem_i=Trem;
%-------------------- BW2 ---------------------------------%     
     [Se, Sp, Sn, Pe, Pp, Pn, CM, Stack]= bw2(e_hat, p_hat, n_hat, Ea, Pa, Na, Sset, M, Stack,Trem);
     
     Tbw2=10;       % Need to find out how much time is spent
    
     Tspent=Tbw1+Tew+Tbw2;
     Trem=Tframe-Tspent;

%------------------------CDTW--------------------------------
    Pe_mat=ones(1,length(Se));         %Every channel gets its own Pe
    Pe_mat=Pe.*Pe_mat;
    Pn_mat=ones(1,length(Sn));          %Every channel gets its own Pn
    Pn_mat=Pn.*Pn_mat;
    Pp_mat=ones(1,length(Sp));         %Every channel gets its own Pp
    Pp_mat=Pp.*Pp_mat;
    [Qe, Qn, Qp, Stack, Prev_set,Ea,Na,Pa,temp_te,temp_tn,temp_tp]=CDTW(CM,Qe,Qn,Qp,Stack,Se,Sp,Sn,Pe_mat,Pp_mat,Pn_mat,Trem, Prev_set,Ea,Na,Pa);
    TP_e(frames)=temp_te/Ne;
    TP_n(frames)=temp_tn/Nn;
    TP_p(frames)=temp_tp/Np;
    
 %-------------------------------Idealised Protocol simulations-----------
 %---We assume BW1 same for both. Further, time is kept same for both-----
    e_hat_i=length(Ea_i);
    if(e_hat_i==0)
        e_hat_i=1;
    end
    n_hat_i=length(Na_i);
    if(n_hat_i==0)
        n_hat_i=1;
    end
    p_hat_i=length(Pa_i);
    if(p_hat_i==0)
        p_hat_i=1;
    end
    %Tspent_i=Tbw1;
    [Se_i, Sp_i, Sn_i, Pe_i, Pp_i, Pn_i, CM_i, Stack_i]= bw2(e_hat_i, p_hat_i, n_hat_i, Ea_i, Pa_i, Na_i, Sset_i, M, Stack_i,Trem_i);
    Tbw2=10;
    Trem_i=Trem_i-Tbw2;
    
 %-----------------------CDTW for idealised protocol-----------------------
    Pe_mat_i=ones(1,length(Se_i));         %Every channel gets its own Pe
    Pe_mat_i=Pe_i.*Pe_mat_i;
    Pn_mat_i=ones(1,length(Sn_i));          %Every channel gets its own Pn
    Pn_mat_i=Pn_i.*Pn_mat_i;
    Pp_mat_i=ones(1,length(Sp_i));         %Every channel gets its own Pp
    Pp_mat_i=Pp_i.*Pp_mat_i;
    [Qe_i, Qn_i, Qp_i, Stack_i, Prev_set_i,Ea_i,Na_i,Pa_i,temp_te_i,temp_tn_i,temp_tp_i]=CDTW(CM_i,Qe_i,Qn_i,Qp_i,Stack_i,Se_i,Sp_i,Sn_i,Pe_mat_i,Pp_mat_i,Pn_mat_i,Trem_i, Prev_set_i,Ea_i,Na_i,Pa_i);
    TP_e_i(frames)=temp_te_i/Ne;
    TP_n_i(frames)=temp_tn_i/Nn;
    TP_p_i(frames)=temp_tp_i/Np;
 end
    index=find(abs(lambda_set-n)<1e-5);
    avgte(RUN,index)=mean(TP_e);
    avgtp(RUN,index)=mean(TP_p);
    avgtn(RUN,index)=mean(TP_n);
    ideal_e(RUN,index)=mean(TP_e_i);
    ideal_p(RUN,index)=mean(TP_p_i);
    ideal_n(RUN,index)=mean(TP_n_i);
    waitbar(n/lambda_set(end),h1,sprintf('current arrival rate= %d%',n));
end
    waitbar(RUN/run_limit,h,sprintf('current run= %d%',RUN));
end

ATe=mean(avgte);
ATp=mean(avgtp);
ATn=mean(avgtn);
ATEi=mean(ideal_e);
ATPi=mean(ideal_p);
ATNi=mean(ideal_n);

total_tp=ATe+ATp+ATn;
ideal_tp=ATEi+ATPi+ATNi;
x_axis=lambda_set;

ATEi=min(ATEi,lambda_set);
ATNi=min(ATNi,lambda_set);
ATPi=min(ATPi,lambda_set);

figure;
plot(x_axis,ATe,'k-o',x_axis,ATn,'k--square',x_axis,ATp,'k:*');
axis tight;
%grid on;
title('Average Throughput v/s \lambda');
xlabel('\lambda');
ylabel('Average Throughput');
legend('Emergency Nodes','Normal Nodes', 'Periodic Nodes');
    
figure;
plot(x_axis,ATe,'k-o',x_axis,ATEi,'k--square',x_axis,x_axis,'k:*');
axis tight;
%grid on;
title('Average Throughput for Emergency Nodes v/s \lambda');
xlabel('\lambda');
ylabel('Average Throughput');
legend('Proposed Protocol','Ideal Protocol','Arrival Rate');

figure;
plot(x_axis,ATn,'k-o',x_axis,ATNi,'k--square',x_axis,x_axis,'k:*');
axis tight;
%grid on;
title('Average Throughput for Normal Nodes v/s \lambda');
xlabel('\lambda');
ylabel('Average Throughput');
legend('Proposed Protocol','Ideal Protocol','Arrival Rate');

figure;
plot(x_axis,ATp,'k-o',x_axis,ATPi,'k--square',x_axis,x_axis,'k:*');
axis tight;
%grid on;
title('Average Throughput for Periodic Nodes v/s \lambda');
xlabel('\lambda');
ylabel('Average Throughput');
legend('Proposed Protocol','Ideal Protocol','Arrival Rate');

figure;
plot(x_axis,ATEi,'k-o',x_axis,ATPi,'k--square',x_axis,ATNi,'k:*');
axis tight;
%grid on;
%title('Average Throughput for Periodic Nodes v/s \lambda');
xlabel('\lambda');
ylabel('Average Throughput');
legend('Ideal emergency','Ideal periodic','Ideal normal');

%save all(1,1,1).mat;





