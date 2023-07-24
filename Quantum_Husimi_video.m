clear all
N=100;

k=sqrt(linspace(1,N-1,N-1));
s=1;
a=diag(k,1);
ac=diag(k,-1);
%gam=0;

Q=sqrt(0.5)*(a+ac);
P=1i*sqrt(0.5)*(ac-a);
q0=5; 
p0=3;
z0=sqrt(0.5)*(q0+p0*1i); %Centre of Gaussian (Amplitude of Coherent State)
omega = 1.0;
gamma = 0.00;
beta = 0.25;
delta = 0;
kappa = 0.0;
epsilon=1.0;

H=(omega-1i*gamma)*(ac*a)+(beta-1i*delta)*(ac*ac*a*a)+(kappa-1i*epsilon)*Q;


%initial state in harmonic oscillator basis:
psi0=zeros(N,1);

psi0(1)=0;
psi0(3)=0;
psi0(4)=1;

psi0=expm(z0*ac-conj(z0)*a)*psi0;


%grid for Husimi plot
D=200;
q=linspace(-10, 10,D);
p=linspace(-10,10,D);

[QH,PH]=meshgrid(q,p);
z=sqrt(0.5)*(QH+i*PH);

%Husimi distribution of initial state:
S=zeros(D,D);
exa=exp(-0.5*abs(z).^2);
for n=1:N
    k=n-1;
    S=S+exa.*conj(z).^k/sqrt(factorial(k))*psi0(n);
end
S=abs(S).^2;

surf(q,p,S,'linestyle','none');
%z = get(h,'ZData');
%set(h,'ZData',z-10);
%caxis([0 1])
view(2)
  axis square
 xlabel('q')
 ylabel('p')

%  
 %time evolution and Husimi functions of time evolution:
 %T=[0.5;2;5;8];
 %T=[pi/4;pi/2;pi];
 T=linspace(0,2*pi,401);
 %T=2*pi;
vidObj=VideoWriter('test3.mp4','MPEG-4')
open(vidObj);
 for ind_t=1:length(T)
     
t=T(ind_t); 
 psi=expm(-1i*H*t)*psi0;

 norm(ind_t)=psi'*psi;
  z_exp(ind_t)=psi'*a*psi/norm(ind_t);
 
%  if mod(ind_t,10)==0
  S=zeros(D,D);
exa=exp(-0.5*abs(z).^2);
for n=1:N
    k=n-1;
    S=S+exa.*conj(z).^k/sqrt(factorial(k))*psi(n);
end
S=abs(S).^2;
 
figure(1)

surf(q,p,S,'linestyle','none');
%z = get(h,'ZData');
%set(h,'ZData',z-10);
%caxis([0 1])
view(2)
  axis square
 xlabel('q')
 ylabel('p')
 title(['t=',num2str(t/pi,'%.2f')])

currFrame = getframe(gcf);
       writeVideo(vidObj,currFrame);
%   hold on 
%   plot(real(z_exp)*sqrt(2),imag(z_exp)*sqrt(2),'w','linewidth',1)
%   pause(0.1)
%   hold off
 end 
 close(vidObj);
% end
% figure
% hold on
%plot(real(z_exp)*sqrt(2),imag(z_exp)*sqrt(2),'m')

 