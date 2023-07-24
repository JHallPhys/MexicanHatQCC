function Classical_Husimi_video
tic
time_final=pi;
%time=[0;0.5;2;8];%
 time=linspace(0,time_final,101);

%grid
gridsize = 10;
q=linspace(-gridsize,gridsize,600);
p=linspace(-gridsize,gridsize,600);
[Q,P]=meshgrid(q,p);
deltaq=q(2)-q(1);
Z=sqrt(0.5)*(Q+i*P);
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
    %Q_prime_ps(ind_p,ind_q,:) = znf(:,3);
    rho0_t(ind_p,ind_q,:)=rho0(zf);
    %rho(ind_p,ind_q,:)=rho0(zf).*Q_prime; %Final solution of the Husimi
    Q_prime_ps(ind_p,ind_q,:) = 1;
    rho(ind_p,ind_q,:)=rho0(zf).*Q_prime; %Final solution of the Husimi
    end
end

%Save and plot Q_prime independently.

vidObj=VideoWriter('Husimi_Medium_0.5Beta_1Ep_PT_symm.mp4','MPEG-4')
open(vidObj);
for ind_t=1:length(t)

figure(1)
surf(q,p,rho(:,:,ind_t),'LineStyle','none')
view(0,90)
axis square
  axis square
 xlabel('q')
 ylabel('p')
 caxis([0,0.03])
 %title(['t=',num2str(t(ind_t)/pi,'%.2f')])
%colorbar

currFrame = getframe(gcf);
       writeVideo(vidObj,currFrame);
% hold on
% plot(q_exp,p_exp,'white');
% xlabel('q')
% ylabel('p')
% pause(0.05)
% hold off
 %   end
end
 close(vidObj);
toc



function dz=deq(t,z)
q=z(1);
p=z(2);
norm=z(3);
omega = 1.0;
gamma = 0.00;
beta = 0.5;
delta = 0.0;
kappa = 0.0;
epsilon = 1;


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
q0=5;
p0=3;
n=0;
z0=sqrt(0.5)*(q0+1i*p0);
r0=abs(z-z0).^(2*n)./factorial(n).*exp(-abs(z-z0).^2);
return

