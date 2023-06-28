clear all
close all
N=120;
D=101;
set_efn='G'; % Invarient Subspace: Gain ('G') Loss ('L')
set_stability='+'; % Stability set: Gain('+'), Loss ('-')
n_efn=6
n0=1
grid_extra=1;
hbar_eff=1;
Xmax=10;
Pmax=10;
set_efn='G'; % Invarient Subspace: Gain ('G') Loss ('L')
set_stability='+'; % Stability set: Gain('+'), Loss ('-')

k=sqrt(linspace(1,N-1,N-1));
s=1;
a=diag(k,1);
ac=diag(k,-1);
%gam=0;
t_final=10
Q=sqrt(0.5)*(a+ac);
P=1i*sqrt(0.5)*(ac-a);
omega = 1.0;
gamma = 0.05;
beta = 0.05;
delta = 1;

H=(-omega-1i*gamma)*(ac*a)+(beta)*(ac*ac*a*a)+(delta)*Q;
U=expm(-1i*H*t_final);
W=U*U';
[psi,En] =eig(W); % psi are the Schur eigenfns and En matrix of eigs
% [psiS,Es]=REig(En,psi,N,set_efn) ;   % Reorder efn/values
Es=diag(En);
[~,ind]=sort(Es,'descend');
Es=Es(ind);
psiS=psi(:,ind);

figure(1)
plot(1:1:N,Es,'k.','Markersize',1)
% return



if ismember(set_efn,'G')

    if ismember(set_stability,'+')
        psi_2=psiS(:,n0:n0+n_efn);
    end
    
      if ismember(set_stability,'-')
         for j =1:n_efn
            psi_2(:,j)=psiS(:,N-n_efn+j);
        end
      end
   
end
    
   
  

q=linspace(-Xmax,Xmax,D);
p=linspace(-Pmax,Pmax,D);

[QH,PH]=meshgrid(q,p);

alpha=sqrt(0.5)*(QH+1i*PH);

Hus_av=zeros(D,D);
Hus=zeros(D,D);
exa=exp(-0.5*abs(alpha).^2);


for j=1:6
    j
% for j=1:20
psi_j=psi_2(:,j);

for n=1:N
    k=n-1;
    Hus=Hus+exa.*conj(alpha).^k/sqrt(factorial(k))*psi_j(n);
%     Hus=Hus+exa.*conj(alpha).^k*coeffs(k+1)*psi_j(n);
end

Hus=abs(Hus).^2;

Hus_av=Hus_av+Hus; % Average 

% Plot each eigenfunction 

figure
clf
imagesc(q,p,Hus);
set(gca,'YDir','normal')
colormap(viridis)
colorbar
% caxis([0 1])
xlabel('q')
ylabel('p') 
Hus(:,:)=0;

end

% Plot average
viridis=viridis();
figure(gcf().Number+1)
clf
imagesc(q*hbar_eff,p*hbar_eff,(Hus_av)) % again scaling of the q p 
set(gca,'YDir','normal')
caxis([0 1])
colormap(viridis)
colorbar
xlabel('q')
ylabel('p')  



sum(sum(Hus_av*abs(q(D)-q(D-1))*abs(p(D)-p(D-1))))  
sum(sum(Hus_av))/(2*pi*n_efn)  

