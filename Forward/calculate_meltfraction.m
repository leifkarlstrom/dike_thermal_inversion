function f=calculate_meltfraction(theta,time,Pm)
%calculate melt fraction based on temperature field 

f=zeros(size(theta));

for i=1:length(time)
    %melt frac in host rock
meltcr=find(theta(Pm.dikepts(end)+1:end,i)>Pm.Tsol&theta(Pm.dikepts(end)+1:end,i)<=Pm.Tliq);
    %melt frac in dike
meltdk=find(theta(1:Pm.dikepts(end),i)>Pm.Tsold&theta(1:Pm.dikepts(end),i)<=Pm.Tliqd);
 if sum(isfinite(meltcr))>0
 f(meltcr,i) = ((theta(meltcr,i)-Pm.Tsol)./(Pm.Tliq-Pm.Tsol)).^Pm.b;
 end
 if sum(isfinite(meltdk))>0
 f(meltdk,i) = ((theta(meltdk,i)-Pm.Tsold)./(Pm.Tliqd-Pm.Tsold)).^Pm.bd;      
 end
 
 if sum(find(f(:,i)>1))>0
     f(f(:,i)>1,i)=1;
 end
end




