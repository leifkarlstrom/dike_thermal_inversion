function TimeTemp = timetemphistory(x,theta,Pts)
%calculate time-temperature history of a point at distances from dike specified by vector Pts
%theta is array with dimensions (distance, time)
%x is spatial vector

TimeTemp = zeros(length(Pts),length(theta(1,:)));

for i=1:length(Pts)
TimeTemp(i,:)=interp1(x,theta,Pts(i));
end

