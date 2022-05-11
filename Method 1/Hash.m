function op = Hash(ip,length)
% Finds Hash Value of a positive integer, Designed to be geometric function
%Ex.- Hash(4)=Hash(100)=0, Hash(15)=Hash(01111)=4.
op_max=ceil(log2(length))+1;
length=op_max;
temp=de2bi(ip,length);
%a=length(temp);
fl=0;
for i=1:length
    if(temp(i)==0)
        fl=1;
        break;
    end
end
if fl==1
    op=i-1;      %Returns position of first zero bit
else
    op=op_max;      %Represents all 1's case
end
    
