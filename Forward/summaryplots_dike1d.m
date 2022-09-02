function summaryplots_dike1d(time,theta,x,Pm)
%requires output from 1D thermal code
%theta = MxN array of temperatures as a function of distance from dike
%center and time
%time, vector of length N

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot time-temperature histories at spatial query points
Pts = linspace(1,100,10);%[2 4 10 40];
%keyboard
for i=1:length(Pts)
TimeTemp(i,:)=interp1(x,theta,Pts(i));
end

figure(1)
for i=1:length(Pts)
    plot(time/Pm.year,TimeTemp(i,:))
    hold on
    xlabel('time in years')
    ylabel('temperature at interogation spatial points (C)')   
    labelT{i}=['x = ' num2str(x(i)) ' m'];
end
set(gca,'xscale','log');
hold off
legend(labelT{1,:})

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%query time points    
int=100;
Tvec=1:int:Pm.TotT;%round(numel(time)/10):numel(time);

%calculate melt fraction associated with temperature
f=calculate_meltfraction(theta,time,Pm);

figure(2)
for i=1:length(Tvec)
    tval=i;
    [~,tval]=min(abs(time-Tvec(i)*Pm.year));
    %subplot(2,1,1)
plot(x,theta(:,tval),'o-')
hold on
xlabel('distance from dike wall in meters')
ylabel('Temperature in C')
title(['Temperature in host rock, curves ploted every ' num2str(int) ' years'])  
% 
%     subplot(2,1,2)
% plot(x,f(:,tval))
% hold on
% xlabel('distance from dike wall in meters')
% ylabel('Melt fraction')
%title(['Melt fraction in host rock, curves ploted every ' num2str(int) ' years'])
end
 hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%calculate max temperature at each spatial point
[MaxTemp, ID]=max(theta,[],2);
MaxTime = time(ID);

figure(3)
scatter(x,MaxTemp,70,MaxTime/Pm.year,'filled','o');colorbar;
xlabel('distance from dike center (m)')
ylabel('Maximum temperature achieved during simulation (C)')
title('Colored by time at which max temperature occurs (years)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%older plotting stuff

%     
% 
% for i=1:length(APts)
% TimeTemp(i,:)=interp1(Pm.x,theta,APts(i));
% FAp(i) = fractionalresetting(TimeTemp(i,:),time,D0a2_A,EaAp,Pm.R);
% FZr(i) = fractionalresetting(TimeTemp(i,:),time,D0a2_Z,EaZr,Pm.R);
% end
% 
% [~,Ti1] = min(abs(time/Pm.year-.007));
% [~,Ti2] = min(abs(time/Pm.year-.9));
% 
% timevals=[time(Ti1:Ti2)/Pm.year logspace(0,4.3010,100)];
% 
% figure(3)
% for i=2:length(APts)%2:20
% 
%     for jj=1:length(timevals)
%         [Tv,Ti] = min(abs(time/Pm.year-timevals(jj)));
%         
%         FracRAp(jj) = fractionalresetting(TimeTemp(i,1:Ti),time(1:Ti),D0a2_A,EaAp,Pm.R);
%         FracRZr(jj) = fractionalresetting(TimeTemp(i,1:Ti),time(1:Ti),D0a2_Z,EaZr,Pm.R);
%     end
%     %i=vec(j);
%     subplot(2,2,2)
%     plot(timevals,FracRAp)
%     set(gca,'xscale','log');
%     xlim([1e-2 2e4])
%     hold on
%     xlabel('time in years')
%     ylabel('partially reset fraction F_d')  
%     subplot(2,2,4)
%     plot(timevals,FracRZr)
%     set(gca,'xscale','log');
%     xlim([1e-2 2e4])
%     hold on
%     xlabel('time in years')
%     ylabel('partially reset fraction F_d') 
%     %label{i}=['x = ' num2str(APts(i)) ' m'];
% end
% set(gca,'xscale','log');
% hold off

   

% % 
% figure(3)
% subplot(2,2,3)
% %hold on
% plot(APts,ApAge_model,'-')
% hold on
% plot(ZPts,ZrAge_model,'-')
% errorbar(Pm.DikeThick + ahe_dist,ahe_age,ahe_err,'o')
% errorbar(Pm.DikeThick + zhe_dist,zhe_age,zhe_err,'o')
% xlabel('distance from contact (m)')
% ylabel('predicted age (Ma)')
% hold off

% clear xM MaxMelt MeltDur
% MeltI=find(sum(f,2)>0);
% for i=1:length(MeltI)
%     IM=MeltI(i);
%     Mdur = find(f(IM,:)>0);
%     durMin=min(Mdur); durMax=max(Mdur);
%     MeltDur(i) = (time(durMax)-time(durMin))/Pm.year;
%     MaxMelt(i) = max(f(IM,:));
%     xM(i)=Pm.x(IM);
% end
% MeltDur(IM+1)=0;
% MaxMelt(IM+1)=0;
% xM(IM+1)=Pm.x(IM+1);
%     
% figure(3)
% %subplot(2,2,4)
% hold off
% yyaxis left
% plot(xM,MaxMelt,'-o');
% ylabel('maximum melt fraction')
% yyaxis right
% plot(xM,MeltDur,'-o');
% ylabel('duration of non-zero melt fraction')
% xlabel('distance from dike center (m)')


% 
% 
% DikeTemp=(1-tanh(10*(time/Pm.year - Pm.tcenter)/Pm.twidth))/2;
% 
% figure(4)
% plot(time/Pm.year,DikeTemp)
% title('dike boundary condition as arctanh')
% ylabel('dike temperature')
% xlim([0 10])
