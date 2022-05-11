function  [aa,bb,cc,dd,ee,ff,gg,hh,ll,mm,nn,oo]=CDTW(CM,Qe,Qn,Qp,Stack,Se,Sp,Sn,Pe_mat,Pp_mat,Pn_mat,Trem,Prev_set,Ea,Na,Pa,time,delay_e,delay_n,delay_p)
% Simulates behavior of CDTW

%Ps_ack=0.8;  %Probability of Successful ack
%Psuc=1;   %Probability of successful tx
Ke=5;      %Upper limit on channel reservation for emergency nodes
Kn=5;      %Upper limit on channel reservation for normal nodes
Kp=5;      %%Upper limit on channel reservation for periodic nodes
Temp_stack=cell(1,1);               %Temparory stack   
% throughput_e=0;
% throughput_n=0;
% throughput_p=0;


%-------------------- Contention Window, Emergency nodes-------------------
N=length(Se);
for i=1:N
    count=0;        % To measure successive empty slots
    Tdtw=0;         % Time spent in data transfer window
    Tcw=0;          % Time spent in contention window
    M=length(CM{Se(i)});
    while ((Tcw<Trem-Tdtw) & (count<3)) % Check for "OR" syntax- checked
        P1=random('bino',1,Pe_mat(i),[1,M]);  %Bernoulli random variable for contention packets
        if sum(P1)>1    % Collision
            Tcw=Tcw+1; 
            count=0;   %to reset count
        elseif sum(P1)==0  % Empty Slot
            Tcw=Tcw+1;
            count=count+1;  % Has to check for consequitivity- done
        else                % Successful Contention
            %Ps= random('bino',1,Ps_ack); % Bernoulli random variable for reservartion packet
            %if Ps==1           %  Check for equal to condition
                index=find(P1,1);   % Finds the first '1' in P1
                node=CM{Se(i)}(index);  %Finds which node is successful
                CM{Se(i)}=setdiff(CM{Se(i)},node);
                tr=min(Ke,length(Qe{node})/2);     %Not sure about syntax- tested ok
                Tdtw=Tdtw+tr;                      %Tr slots reserved
                Stack{Se(i)}=[Stack{Se(i)} [node,tr]]; %Stack saved for Data transfer window
                Tcw=Tcw+1;
                count=0;      %to reset count
                Pe_mat(i)=1/((1/Pe_mat(i))-1);           %Update Pe after a successful contention
                Pe_mat(i)=min(Pe_mat(i),0.8);
            %else
                %count=count+1;  % Loss of contention packet equivalent to empty slot for BS
                %Tcw=Tcw+1;           % Loss of Contention Packet
            %end
        end
        M=length(CM{Se(i)});        % No of users in "i'th" channel
    end
    time_fr=time+Tcw;
    Tsilent=Trem-(Tcw+Tdtw);
    time_fr=time_fr+Tsilent;
    % Write code for Data transfer window now
    
    %---------------------------------Data Transfer Window, Emergency node---------------
    S_l=length(Stack{Se(i)});       %Size of the Stack
       
    while S_l~=0
        Temp=Stack{Se(i)}(S_l-1:S_l);          %Stack popped
        Stack{Se(i)}(S_l-1:S_l)=[];            %Last two elements removed
        node=Temp(1);
        slots=Temp(2);
        %Ps= random('bino',1,Psuc,[1,slots]);   %Array of successful tx
        %total=slots;                   %Total successful tx
        %throughput_e=throughput_e+total;          %No of successful tx
        for tr_time=1:slots
            arrival=Qe{node}(2);
            time_fr=time_fr+1;
            delay_e{node}=[delay_e{node} [time_fr-arrival]];
            Qe{node}(1:2)=[];
        end
        
        %Qe(node)=Qe(node)-total;        % Queue of the node updated
        if(isempty(Qe{node}))
            Ea=setdiff(Ea,node);           %Removing the node from Ea if the queue becomes zero
        end
        S_l=length(Stack{Se(i)});
    end
end

%-- Contention Window, normal nodes------------------------

N=length(Sn);
for i=1:N
    count=0;        % To measure successive empty slots
    Tdtw=0;         % Time spent in data transfer window
    Tcw=0;          % Time spent in contention window
    M=length(CM{Sn(i)});
    while ((Tcw<Trem-Tdtw) & (count<3)) % Check for "OR" syntax- checked
        P1=random('bino',1,Pn_mat(i),[1,M]);  %Bernoulli random variable for contention packets
        if sum(P1)>1    % Collision
            Tcw=Tcw+1; 
            count=0;   %to reset count
        elseif sum(P1)==0  % Empty Slot
            Tcw=Tcw+1;
            count=count+1;  % Has to check for consequitivity- done
        else                % Successful Contention
            %Ps= random('bino',1,Ps_ack); % Bernoulli random variable for reservartion packet
            %if Ps==1           %  Check for equal to condition
                index=find(P1,1);   % Finds the first '1' in P1
                node=CM{Sn(i)}(index);  %Finds which node is successful
                CM{Sn(i)}=setdiff(CM{Sn(i)},node);
                tr=min(Kn,length(Qn{node})/2);     %Not sure about syntax- tested ok
                Tdtw=Tdtw+tr;                      %Tr slots reserved
                Stack{Sn(i)}=[Stack{Sn(i)} [node,tr]]; %Stack saved for Data transfer window
                Tcw=Tcw+1;
                count=0;      %to reset count
                Pn_mat(i)=1/((1/Pn_mat(i))-1);           %Update Pn after a successful contention
                Pn_mat(i)=min(Pn_mat(i),0.8);
            %else
                %count=count+1;  % Loss of contention packet equivalent to empty slot for BS
                %Tcw=Tcw+1;           % Loss of Contention Packet
            %end
        end
        M=length(CM{Sn(i)});        % No of users in "i'th" channel
    end
    % Write code for Data transfer window now
    time_fr=time+Tcw;
    Tsilent=Trem-(Tcw+Tdtw);
    time_fr=time_fr+Tsilent;
    %---------------------------------Data Transfer Window, normal node---------------
    S_l=length(Stack{Sn(i)});       %Size of the Stack
       
    while S_l~=0
        Temp=Stack{Sn(i)}(S_l-1:S_l);          %Stack popped
        Stack{Sn(i)}(S_l-1:S_l)=[];            %Last two elements removed
        node=Temp(1);
        slots=Temp(2);
        %Ps=random('bino',1,Psuc,[1,slots]);    %Array to represent transmitts     
        %total=slots;                 %total successful tx
        %throughput_n=throughput_n+total;
        
        for tr_time=1:slots
            arrival=Qn{node}(2);
            time_fr=time_fr+1;
            delay_n{node}=[delay_n{node} [time_fr-arrival]];
            Qn{node}(1:2)=[];
        end
        %Qn(node)=Qn(node)-total;        % Queue of the node updated
        if(isempty(Qe{node}))
            Na=setdiff(Na,node);           %Removing the node from Na if the queue becomes zero
        end
        S_l=length(Stack{Sn(i)});
    end
end

%------------Contention Window, periodic nodes---------------------------------

N=length(Sp);
for i=1:N
    count=0;        % To measure successive empty slots
    Tdtw=0;         % Time spent in data transfer window
    Tcw=0;          % Time spent in contention window
    S_l=length(Stack{Sp(i)});       %Size of the Stack
    if S_l~=0                   %Check for the nodes who have reserved access in previous frames
        prev_slots=S_l/2;        %Don't need toactuaaly pop the whole stack as length of stack is readily known.    
        Tdtw=Tdtw+prev_slots;
    end
    M=length(CM{Sp(i)});  
    while ((Tcw<Trem-Tdtw) & (count<3)) % Check for "OR" syntax- checked
        P1=random('bino',1,Pp_mat(i),[1,M]);  %Bernoulli random variable for contention packets
        if sum(P1)>1    % Collision
            Tcw=Tcw+1; 
            count=0;   %to reset count
        elseif sum(P1)==0  % Empty Slot
            Tcw=Tcw+1;
            count=count+1;  % Has to check for consequitivity- done
        else                % Successful Contention
            %Ps= random('bino',1,Ps_ack); % Bernoulli random variable for reservartion packet
            %if Ps==1           %  Check for equal to condition
                index=find(P1,1);   % Finds the first '1' in P1
                node=CM{Sp(i)}(index);  %Finds which node is successful
                CM{Sp(i)}=setdiff(CM{Sp(i)},node);
                fr=min(Kp,length(Qp{node,1})/2);     %Not sure about syntax- tested ok
                Tdtw=Tdtw+1;                      %tr slots reserved
                Stack{Sp(i)}=[Stack{Sp(i)} [node,fr]]; %Stack saved for Data transfer window
                Tcw=Tcw+1;
                count=0;      %to reset count
                Pp_mat(i)=1/((1/Pp_mat(i))-1);           %Update Pp after a successful contention
                Pp_mat(i)=min(Pp_mat(i),0.8);
            %else
                %count=count+1;  % Loss of contention packet equivalent to empty slot for BS
                %Tcw=Tcw+1;           % Loss of Contention Packet
            %end
        end
        M=length(CM{Sp(i)});        % No of users in "i'th" channel
    end
    % Write code for Data transfer window now
    time_fr=time+Tcw;
    Tsilent=Trem-(Tcw+Tdtw);
    time_fr=time_fr+Tsilent;
    %---------------------------------Data Transfer Window, Periodic node---------------
    S_l=length(Stack{Sp(i)});       %Size of the Stack
       
    while S_l~=0
        Temp=Stack{Sp(i)}(S_l-1:S_l);          %Stack popped
        Stack{Sp(i)}(S_l-1:S_l)=[];            %Last two elements removed
        node=Temp(1);
        frames=Temp(2);
%         slots=1;                               %Since Periodic nodes transmit data in 1 slot only
%         for j=1:slots
             %Ps= random('bino',1,Psuc);     
%              Ps=1;
%              if Ps==1                      %Successful TX
                 %throughput_p=throughput_p+1;          %No of successful tx, periodic
                 %Qp(node,1)=Qp(node,1)-1;        % Queue of the node updated
                 
                 arrival=Qp{node,1}(2);
                 time_fr=time_fr+1;
                 delay_p{node}=[delay_p{node} [time_fr-arrival]];
                 Qp{node,1}(1:2)=[];
            
                 frames=frames-1;
                 if(frames>0)            %Reserved access in next frame also
                     Temp_stack{1}=[Temp_stack{1} [node,frames]];
                     Prev_set=union(Prev_set,node);    %Node added to Prev_set   
                     Pa=setdiff(Pa,node);           %Removed the set from Pa so it doesn't participate in further CWs
                 else                        %Access finished
                    Prev_set=setdiff(Prev_set,node);     %Prev_set updated 
                 end
                 if(isempty(Qp{node,1}))         %Queue is empty
                     if(isempty(Qp{node,2}))
                         Pa=setdiff(Pa,node);           %Removed the set from Pa
                     end
                     Qp{node,1}=Qp{node,2};
                     Qp{node,2}=[];                    %Queue of the node updated
                 end
%              else                %Unsuccessful Tx
%                  frames=frames-1;
%                  if(frames>0)       %Access reserved in next frame also
%                     Temp_stack{1}=[Temp_stack{1} [node,frames]];
%                     Prev_set=union(Prev_set,node);    %Node added to Prev_set   
%                     Pa=setdiff(Pa,node);           %Removed the set from Pa so it doesn't participate in further CWs
%                  else            %Access finished
%                      Prev_set=setdiff(Prev_set,node);     %Prev_set updated 
%                  end             
%              end
%         end
        S_l=length(Stack{Sp(i)});
    end
    S_lt=length(Temp_stack{1});                 % Restoring Stack from temp_stack
    while S_lt~=0
        Temp1=Temp_stack{1}(S_lt-1:S_lt);          %Stack popped
        Temp_stack{1}(S_lt-1:S_lt)=[];            %Last two elements removed
        Stack{Sp(i)}=[Stack{Sp(i)} Temp1];       %Stack updated
        S_lt=length(Temp_stack{1});
    end
end
time=time+Trem;

aa=Qe;
bb=Qn;
cc=Qp;
dd=Stack;
ee=Prev_set;
ff=Ea;
gg=Na;
hh=Pa;
%ii=throughput_e;
%jj=throughput_n;
%kk=throughput_p;
mm=delay_e;
ll=time;
nn=delay_n;
oo=delay_p;
                
            
            
    