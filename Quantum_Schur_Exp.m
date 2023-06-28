clear all
close all
N=51;
D=51;
set_efn='G'; % Invarient Subspace: Gain ('G') Loss ('L')
set_stability='+'; % Stability set: Gain('+'), Loss ('-')
n_efn=2
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

Q=sqrt(0.5)*(a+ac);
P=1i*sqrt(0.5)*(ac-a);
omega = 1.0;
gamma = 0.05;
beta = 0.05;
delta = 1;

H=(-omega-1i*gamma)*(ac*a)+(beta)*(ac*ac*a*a)+(delta)*Q;
[psi,En] = schur(H); % psi are the Schur eigenfns and En matrix of eigs
[psiS,Es]=REig(En,psi,N,set_efn) ;   % Reorder efn/values
Es=diag(Es);
if ismember(set_efn,'G')

    if ismember(set_stability,'+')
        psi_2=psiS(:,1:n_efn);
    end
    
      if ismember(set_stability,'-')
         for j =1:n_efn
            psi_2(:,j)=psiS(:,N-n_efn+j);
        end
      end
   
end
    

psi0=psiS(:,1);
Xop=Q; % Position operator
EX=psi0'*Xop*psi0
EXS=psi0'*Xop^2*psi0
sigmax=real(sqrt(EXS-EX^2)) % sigma_x

Pop=P; % Momentum operator
EP=psi0'*Pop*psi0;
EPS=psi0'*Pop^2*psi0;
sigmap=real(sqrt(EPS-EP^2)) % Sigma p 

sigmax*sigmap
return


% figure(2)
% hold on
% plot(1:1:N,imag(Es),'b.')
% return
% axis([-250 50 -80 0])
% xlabel('Re(\lambda)')
% ylabel('Im(\lambda)')

% % return
for j = 1:100
    Qexp(j)=psiS(:,j)'*Q*psiS(:,j);
    Pexp(j)=psiS(:,j)'*P*psiS(:,j);
end


figure(1)  
clf
l1=plot(real(Qexp),real(Pexp),'g.','Markersize',5)
   
  



%==========================================================================
%==========================================================================
%==========================================================================
% Its a sin
%==========================================================================
%==========================================================================

system_flag='MKS';
time_flag='fwd';
mu=0.001 % Parameter that measures distance to fp at t=0
N0=1; % Number of points
t_final=80;
t_final_bwd=30

qmin=-8
qmax=8
pmin=-8
pmax=8
% qmin=-20
% qmax=20
% pmin=-20
% pmax=20


nfp=2
%==========================================================================
% Construct Fixed points
%==========================================================================

[q_fp,p_fp]=get_fixed_points_mexican_hat(omega,beta,delta,gamma,system_flag); 



%==========================================================================
% Construct the Jacobian
%==========================================================================


[J_out,Evec,Eval]=get_jacobian_mexican_hat(q_fp(nfp),p_fp(nfp),omega,beta,delta,gamma,system_flag); 

q_fp=q_fp(nfp);
p_fp=p_fp(nfp);

figure(1)
hold on
l2=plot(q_fp,p_fp,'k.','Markersize',5)
xlabel('q')
ylabel('p')
axis([qmin qmax pmin pmax])




%==========================================================================
% Set up the initial conditions
%==========================================================================

Eval=diag(Eval);
[~,ind_growing]=max(Eval); % Index of growing eignvalues

u1=Evec(:,1);
u1=u1./norm(u1); % Normnalised eigenvector in growing direction

u2=Evec(:,2);
u2=u2./norm(u2); % Normnalised eigenvector in growing direction



q_0=linspace(0,1,N0);
p_0=linspace(0,1,N0);


% 
% q_stable=q_fp+q_0*mu*u1(1);
% p_stable=p_fp+p_0*mu*u1(2);

%==========================================================================
%==========================================================================
% This is bad and needs replaced
% Know shouldnt be doing thiz
%==========================================================================
%==========================================================================


% One side

q_stable=q_fp-q_0*mu*u2(1);
p_stable=p_fp-p_0*mu*u2(2);



for j=1:length(q_stable)
    j
    % Counter for the number we are on

        % 1) Initial conditions in phase space       
        z(1)=q_stable(j);
        z(2)=p_stable(j);
        z(3)=1;
        
          % 2) Propagate Norm
        
   [t,z_ps] = ode89(@(t,z) dhat(t,z,omega,beta,delta,gamma,system_flag,time_flag) ,[0 t_final],z); % Integrate
        
    figure(1)
    hold on
    plot(z_ps(:,1),z_ps(:,2),'r-','Markersize',1)
    xlabel('q')
    ylabel('p')
    axis([qmin qmax pmin pmax])

    pause(1)
       
    
end


% Other side

q_stable=q_fp+q_0*mu*u2(1);
p_stable=p_fp+p_0*mu*u2(2);


%==========================================================================
% Okay work the manifold
%==========================================================================


for j=1:length(q_stable)
    j
    % Counter for the number we are on

        % 1) Initial conditions in phase space       
        z(1)=q_stable(j);
        z(2)=p_stable(j);
        z(3)=1;
        
          % 2) Propagate Norm
        
       [t,z_ps] = ode89(@(t,z) dhat(t,z,omega,beta,delta,gamma,system_flag,time_flag) ,[0 t_final],z); % Integrate
        
    figure(1)
    hold on
    l3=plot(z_ps(:,1),z_ps(:,2),'r-','Markersize',1)
    xlabel('q')
    ylabel('p')
    axis([qmin qmax pmin pmax])

    pause(1)
       
    
end

% return

%==========================================================================
%==========================================================================
% Untable manifold
%==========================================================================
%==========================================================================

time_flag='bwd';
t_final=t_final_bwd
% One side

q_stable=q_fp+q_0*mu*u1(1);
p_stable=p_fp+p_0*mu*u1(2);



for j=1:length(q_stable)
    j
    % Counter for the number we are on

        % 1) Initial conditions in phase space       
        z(1)=q_stable(j);
        z(2)=p_stable(j);
        z(3)=1;
        
          % 2) Propagate Norm
        
      [t,z_ps] = ode89(@(t,z) dhat(t,z,omega,beta,delta,gamma,system_flag,time_flag) ,[0 t_final],z); % Integrate
        
        
    figure(1)
    hold on
    plot(z_ps(:,1),z_ps(:,2),'b-','Markersize',1)
    xlabel('q')
    ylabel('p')
    axis([qmin qmax pmin pmax])

    pause(1)
       
    
end


% Other side

q_stable=q_fp-q_0*mu*u1(1);
p_stable=p_fp-p_0*mu*u1(2);



for j=1:length(q_stable)
    j
    % Counter for the number we are on

        % 1) Initial conditions in phase space       
        z(1)=q_stable(j);
        z(2)=p_stable(j);
        z(3)=1;
        
          % 2) Propagate Norm
        
   [t,z_ps] = ode89(@(t,z) dhat(t,z,omega,beta,delta,gamma,system_flag,time_flag) ,[0 t_final],z); % Integrate
        
        
    figure(1)
    hold on
    l4=plot(z_ps(:,1),z_ps(:,2),'b-','Markersize',1)
    xlabel('q')
    ylabel('p')
    legend([l1 l2 l3 l4],{'Eigenstate Expectation','Fixed point','Stable manifold','Unstable manifold'})
    axis([qmin qmax pmin pmax])
    box on
    pause(1)
       
    
end
