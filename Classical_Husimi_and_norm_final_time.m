function Classical_Husimi_and_norm_final_time
time_final=8;
tic

%grid
gridsize = 10;
q=linspace(-gridsize,gridsize,200);
p=linspace(-gridsize,gridsize,200);
[Q,P]=meshgrid(q,p);
dq = abs(q(2) - q(1));
dp = abs(p(2) - p(1));

n0=1;
for ind_q=1:length(q) %index
q_ic=q(ind_q); %initial condition for this ODE run
    for ind_p=1:length(p)
    p_ic=p(ind_p);%initial condition for this ODE run
    [t,znf]=ode45(@(t,z) deq(t,z),[0 time_final],[q_ic,p_ic,n0]); %output evolved phase-space points at a time as array
    zf(ind_p,ind_q)=sqrt(0.5)*(znf(end,1)+1i*znf(end,2)); %constructing single vector of phase-space points
  W(ind_p,ind_q)=znf(end,3);
    % Q_prime=znf(:,3); %The solution to the dQ/dt OD 
    end
end

% zf = zf./sqrt(2);
% 
% q_val_1 = 1.0563;
% p_val_1 = -0.0559;
% 
% plot(zf(:,1),zf(:,2));
% hold on
% plot(q_val_1,q_val_2)


    rho_f=rho0(zf).*W; %Final solution of the Husimi
   husimi_norm= sum(sum(rho_f))*dp*dq;
rho_f=rho_f/husimi_norm;
%sum(sum(rho_f))*dp*dq


% figure(1)
% surf(q,p,rho_f,'LineStyle','none')
% %hold on
% %plot(5*sin(theta), 5*cos(theta),'red');
% xlabel('q')
% ylabel('p')
% view(0,90)
% axis square
% colorbar


figure(5)
surf(q,p,real(log(real(W))),'LineStyle','none')
%hold on
%plot(5*sin(theta), 5*cos(theta),'red');
xlabel('q')
ylabel('p')
view(0,90)
axis square
colorbar
%end
title('"Norm Husimi Only"')
caxis([-80 0])

figure(6)
surf(q,p,rho0(zf),'LineStyle','none')
%hold on
%plot(5*sin(theta), 5*cos(theta),'red');
xlabel('q')
ylabel('p')
view(0,90)
axis square
colorbar
title('"Q prime Husimi Only"')
%end
figure(7)
surf(q,p,rho_f,'LineStyle','none')
%hold on
%plot(5*sin(theta), 5*cos(theta),'red');
xlabel('q')
ylabel('p')
view(0,90)
axis square
colorbar
caxis([0 0.08])
toc;

% vidObj=VideoWriter('Norm_CHO_classical.mp4','MPEG-4')
% open(vidObj);
% 
% for ind_t=1:length(t)
% 
% figure(1)
% surf(q,p,log(real(W(:,:,ind_t))),'LineStyle','none')
% view(0,90)
% axis square
%   axis square
%  xlabel('q')
%  ylabel('p')
%  title(['t=',num2str(t(ind_t)/pi,'%.2f')])
% 
% currFrame = getframe(gcf);
%        writeVideo(vidObj,currFrame);
% end
%  close(vidObj);


function dz=deq(t,z)
q=z(1);
p=z(2);
norm=z(3);
omega = -1.0;
gamma = 0.05;
beta = 0.05;
delta = 0.0;
kappa = 1.0;
epsilon = 0.0;


%% Miburn, (w-igamma)HO + (beta -idelta)HO^2 + (Kappa -i epsilon)q
%H = 0.5*(omega-1i*gamma)*(p^2+q^2) +0.25*(beta-1i*delta)*(p^2+q^2)^2 ...
%    + (kappa - 1i*epsilon)*q;
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
