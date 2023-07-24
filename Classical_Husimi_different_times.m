function Classical_Husimi_different_times
clear all
%time=[0;1;2;5];%
plot_time=[0;0.5;2;8];
%plot_time=8;
time=linspace(0,max(plot_time),101);
 

%grid
gridsize = 10;
q=linspace(-gridsize,gridsize,200); %240
p=linspace(-gridsize,gridsize,200);
[Q,P]=meshgrid(q,p);
deltaq=q(2)-q(1);
Z=sqrt(0.5)*(Q+1i*P);
n0=1;
rho=NaN*ones(length(p),length(q),length(time));
rho0_t=NaN*ones(length(p),length(q),length(time));
%Q_prime=NaN*zeros(length(time));
Q_prime_ps=NaN*ones(length(p),length(q),length(time));


q_exp=sum(sum(Q.*rho0(Z)))*deltaq^2/2/pi;
p_exp=sum(sum(P.*rho0(Z)))*deltaq^2/2/pi;

for ind_q=1:length(q) %index
q_ic=q(ind_q); %initial condition for this ODE run
    for ind_p=1:length(p)
    p_ic=p(ind_p);%initial condition for this ODE run
    [t,znf]=ode45(@(t,z) deq(t,z),time,[q_ic,p_ic,n0]); %output evolved phase-space points at a time as array
    zf=sqrt(0.5)*(znf(:,1)+1i*znf(:,2)); %constructing single vector of phase-space points
    Q_prime=znf(:,3); %The solution to the dQ/dt ODE
    Q_prime_ps(ind_p,ind_q,:) = znf(:,3);
    rho0_t(ind_p,ind_q,:)=rho0(zf);
    rho(ind_p,ind_q,:)=rho0(zf).*Q_prime; %Final solution of the Husimi
    end
end

%Save and plot Q_prime independently.

for ind_t=1:length(t)
% 
 q_exp(ind_t)=sum(sum((Q.*rho(:,:,ind_t))))*deltaq^2/2/pi;
 p_exp(ind_t)=sum(sum((P.*rho(:,:,ind_t))))*deltaq^2/2/pi;
 norm(ind_t)=sum(sum(rho(:,:,ind_t)))*deltaq^2/2/pi;
 q_exp(ind_t)=q_exp(ind_t)/norm(ind_t);
 p_exp(ind_t)=p_exp(ind_t)/norm(ind_t);

 dq = abs(q(2) - q(1));
 dp = abs(p(2) - p(1));
 husimi_norm= sum(sum(rho))*dp*dq;
 rho=rho./husimi_norm;

%% Plots Remove if loop and copies if want a video
   if abs(time(ind_t)) == 0
f= figure(1)
%exportsetupdlg(f)
surf(q,p,rho(:,:,ind_t)-1,'LineStyle','none')
view(0,90)
axis square
%colorbar
%title(['t=',num2str(t(ind_t))])
%title('t = \pi/2')
hold on
plot(q_exp,p_exp,'white');
xlabel('q')
ylabel('p')
pause(1)
hold off
   end
   if abs(time(ind_t) - plot_time(2)) <max(plot_time)/101
f = figure(2)
%exportsetupdlg(f)
surf(q,p,rho(:,:,ind_t),'LineStyle','none')
view(0,90)
axis square
%colorbar
%title(['t=',num2str(t(ind_t))])
%title('t = \pi')
hold on
%plot(q_exp,p_exp,'white');
xlabel('q')
ylabel('p')
pause(1)
hold off
     end
     if abs(time(ind_t) - plot_time(3)) <max(plot_time)/101
f = figure(3)
%exportsetupdlg(f)
surf(q,p,rho(:,:,ind_t),'LineStyle','none')
view(0,90)
axis square
%colorbar
%title(['t=',num2str(t(ind_t))])
%title('t = \pi')
hold on
%plot(q_exp,p_exp,'white');
xlabel('q')
ylabel('p')
pause(1)
hold off
     end
 if time(ind_t) == plot_time(4)
f = figure(4)
%exportsetupdlg(f)
surf(q,p,rho(:,:,ind_t),'LineStyle','none')
%surf(q,p,rho(:,:,ind_t)-1,'LineStyle','none')
view(0,90)
axis square
%colorbar
%title(['t=',num2str(t(ind_t))])
%title('t = 2\pi')
hold on
%plot(q_exp,p_exp,'white');
xlabel('q')
ylabel('p')
xlim([-10 10])
ylim([-10 10])
pause(1)
hold off
end
end
% title('"Full Husimi"')
% 
% 
% %hold off
% figure(9)
% surf(q,p,rho0_t(:,:,length(t)),'LineStyle','none')
% %hold on
% %plot(5*sin(theta), 5*cos(theta),'red');
% xlabel('q')
% ylabel('p')
% view(0,90)
% axis square
% colorbar
% title('"Evolved Initial Husimi Only"')
% 
% %for ind_t=1:length(t)
%      ind_t=length(t);
% %hold off
% figure(2)
% surf(q,p,Q_prime_ps(:,:,ind_t),'LineStyle','none')
% %hold on
% %plot(5*sin(theta), 5*cos(theta),'red');
% xlabel('q')
% ylabel('p')
% view(0,90)
% axis square
% colorbar
% %end
% %title('"Norm Husimi Only"')




function dz=deq(t,z)
q=z(1);
p=z(2);
norm=z(3);
omega = -1.0;%
gamma = 0.05; %
beta = 0.05;%
delta = 0.0;
kappa = 1.0;%
epsilon = 0.0;


%% Miburn, (w-igamma)HO + (beta -idelta)HO^2 + (Kappa -i epsilon)q
%H = 0.5*(omega-1i*gamma + beta)*(p^2+q^2) +0.25*(beta-1i*delta)*(p^2+q^2)^2 ...
%    + (kappa - 1i*epsilon)*q; for the PT symm and Mexican hat systems
 G = -0.5*gamma*(p.^2+q.^2) -0.25*delta*(p.^2+q.^2).^2 - epsilon*q;
 Hq= (omega+beta)*q + beta*(p.^2+q.^2).*q + kappa;
 Hp= (omega+beta)*p + beta*(p.^2+q.^2).*p;
 Gq=-gamma*q - delta*(p.^2+q.^2).*q - epsilon;
 Gp=-gamma*p - delta*(p.^2+q.^2).*p;


%% For MoC trajectories
dz=[-Hp+Gq;
    Hq+Gp;
    2*G*norm];

return


function r0=rho0(z);
q0=-3;
p0=5;
n=2;
z0=sqrt(0.5)*(q0+1i*p0);
r0=abs(z-z0).^(2*n)./factorial(n).*exp(-abs(z-z0).^2);
return

