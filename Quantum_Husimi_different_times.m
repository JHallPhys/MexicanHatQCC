clear all
N=200;

k=sqrt(linspace(1,N-1,N-1));
s=1;
a=diag(k,1);
ac=diag(k,-1);
%gam=0;

Q=sqrt(0.5)*(a+ac);
P=1i*sqrt(0.5)*(ac-a);
q0=-3; 
p0=5;
z0=sqrt(0.5)*(q0+p0*1i); %Centre of Gaussian (Amplitude of Coherent State)
omega = -1.0;%
gamma = 0.05; %
beta = 0.05;%
delta = 0.0;
kappa = 1.0;%
epsilon = 0.0;


H=(omega-1i*gamma)*(ac*a)+(beta-1i*delta)*(ac*ac*a*a)+(kappa-1i*epsilon)*Q;
%H=(omega-1i*gamma + beta)*(ac*a)+(beta-1i*delta)*(ac*ac*a*a)+(kappa-1i*epsilon)*(ac + a)/sqrt(2);

%initial state in harmonic oscillator basis:
psi0=zeros(N,1);

psi0(1)=0;
%psi0(1) = 1;
psi0(3)=1; %third state is psi(4)

psi0=expm(z0*ac-conj(z0)*a)*psi0;


%grid for Husimi plot
D=450;
q=linspace(-10, 10,D);
p=linspace(-10,10,D);
dq = abs(q(2)-q(1));
dp = abs(p(2)-p(1));

[QH,PH]=meshgrid(q,p);
z=sqrt(0.5)*(QH+1i*PH);

%Husimi distribution of initial state:
S=zeros(D,D);
exa=exp(-0.5*abs(z).^2);
for n=1:N
    k=n-1;
    S=S+exa.*conj(z).^k/sqrt(factorial(k))*psi0(n);
end
S=abs(S).^2;

husimi_norm =sum(sum(S))*dp*dq; %equals 2pi, ideal for Husimi.
S = S/husimi_norm;
sum(sum(S))*dp*dq; %equals 1, shows renormalised Husimi
 
 %time evolution and Husimi functions of time evolution:
 T=[0;0.5;2;8];
 %T=[0;pi/10;pi/4;pi];
 T_1=linspace(0,max(T),201);
 %T=2*pi;

for ind_t=1:length(T_1)
t=T_1(ind_t); 
 psi=expm(-1i*H*t)*psi0;

 norm(ind_t)=psi'*psi;
  z_exp(ind_t)=psi'*a*psi/norm(ind_t);
end

 for ind_t=1:length(T)
     t=T(ind_t); 
    psi=expm(-1i*H*t)*psi0;
%  if mod(ind_t,10)==0
 state_norm=psi'*psi;
 psi = psi/sqrt(state_norm);
  S=zeros(D,D);
exa=exp(-0.5*abs(z).^2);
for n=1:N
    k=n-1;
    S=S+exa.*conj(z).^k/sqrt(factorial(k))*psi(n);
end
S=abs(S).^2;
S = S/(2*pi);
sum(sum(S))*dp*dq; %sanity check, equals 1

f= figure(ind_t)
%exportsetupdlg(f)
%figure(1)

surf(q,p,S,'linestyle','none');
%z = get(h,'ZData');
%set(h,'ZData',z-10);
%caxis([0 1])
view(2)
  axis square
 xlabel('q')
 ylabel('p')
 %title(['t=',num2str(t)])
%  if ind_t ==2
%     title('t = \pi/2')
%         elseif ind_t ==3
%         title('t = \pi')
%             elseif ind_t == 4
%             title('t = 2\pi')
%  end
  %hold on 
  %for ind_t_1=1:length(T_1) %plot
  %    if T_1(ind_t_1) <= t
  %        z_exp_1(ind_t_1) = z_exp(ind_t_1);
  %    end
  %end
  %plot(real(z_exp_1)*sqrt(2),imag(z_exp_1)*sqrt(2),'w','linewidth',1)
  %pause(0.05)
  %hold off
 end 
 hold off
% end
% figure
%hold on
%plot(real(z_exp)*sqrt(2),imag(z_exp)*sqrt(2),'m')

 