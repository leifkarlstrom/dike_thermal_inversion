function Age_model = thermochron_calc(data_dist,theta,x,time,ThermParam,Ea,Tr,Ti,Tm,Pm,flag)
%keyboard
TcPts = Pm.DikeThick + data_dist*Pm.dip; %evaluation points for Zr

%calculate time-temperature histories at query points
TimeTemp = timetemphistory(x,theta,TcPts);

F=zeros(length(TcPts),1);
for i=1:length(TcPts)
%fractional resetting of Zr and Ap
    if flag==1 %spherical grain
        D0a2=ThermParam;        
        F(i) = fractionalresetting_sphere(TimeTemp(i,:),time,D0a2,Ea,Pm.R);
    elseif flag==2 %cylindrical grain
        D0a2=ThermParam;
        F(i) = fractionalresetting_cylinder(TimeTemp(i,:),time,D0a2,Ea,Pm.R);

    elseif flag==3 %apatite fission track
        rmr0=ThermParam;
        F(i) = fractionalresetting_FT(TimeTemp(i,:),time,Tm,rmr0,Pm.year);        
% %zft
%     else flag==4;
%         F(i) = fractionalresetting_FT(TimeTemp(i,:),time,Tm);
    end
end
%predicted age based on fractional resetting
Age_model = Tr + (Ti-Tr).*(1-F);

if flag==3
    %apply standard length correction from Ketcham to age
    Age_model = Age_model/0.893;
end

