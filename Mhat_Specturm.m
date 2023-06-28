clear all
close all
N=151;
D=51;
set_efn='G'; % Invarient Subspace: Gain ('G') Loss ('L')
set_stability='+'; % Stability set: Gain('+'), Loss ('-')
n_efn=1
grid_extra=1;
hbar_eff=1;
Xmax=10;
Pmax=10;


k=sqrt(linspace(1,N-1,N-1));
s=1;
a=diag(k,1);
ac=diag(k,-1);
%gam=0;

Q=sqrt(0.5)*(a+ac);
P=1i*sqrt(0.5)*(ac-a);
omega = 1.0;
gamma = 0.05;
beta = 0.05;
delta = 0.1;
 'LOL'

H=(-omega-1i*gamma)*(ac*a)+(beta)*(ac*ac*a*a)+(delta)*Q;
[psi,En] = schur(H); % psi are the Schur eigenfns and En matrix of eigs
[psiS,Es]=REig(En,psi,N,set_efn) ;   % Reorder efn/values

Es=diag(Es);

figure(1)
hold on
plot(real(Es),imag(Es),'b.','Markersize',10)
xlabel('Re(\lambda)')
ylabel('Im(\lambda)')

figure(2)
hold on
plot(1:1:N,imag(Es),'b.','Markersize',10)
xlabel('Eigenvalue Index')
ylabel('Im(\lambda)')
return
figure(3)
hold on
plot(real(Es(1:4)),imag(Es(1:4)),'ro','Markersize',2)
plot(real(Es(5:19)),imag(Es(5:19)),'gx','Markersize',2)
plot(real(Es(20:25)),imag(Es(20:25)),'b.','Markersize',2)
xlabel('Re(\lambda)')
ylabel('Im(\lambda)')
legend('Branch 1','Branch 2','Branch 3','Location','northwest')

% Level spacng

ds_imag=zeros(N-1,1);
for j=1:N-1
    ds_imag(j)=abs(imag(Es(j+1))-imag(Es(j)));
end

figure
plot(1:1:N-1,ds_imag,'k.')
xlabel('n')
ylabel('spacing')

return

