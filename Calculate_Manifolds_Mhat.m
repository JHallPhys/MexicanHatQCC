clear all 
close all
%==========================================================================
% System parameters
%==========================================================================

omega=1;%
gamma=0.05;
beta=0.05;
delta=1;
system_flag='MKS';
time_flag='fwd';
mu=0.001 % Parameter that measures distance to fp at t=0
t_final=80;
t_final_bwd=30 % If this is too high unstable manifold will move out too far
qmin=-8
qmax=8
pmin=-8
pmax=8
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
l1=plot(q_fp,p_fp,'k.','Markersize',5)
xlabel('q')
ylabel('p')
axis([qmin qmax pmin pmax])

%==========================================================================
% Set up the initial conditions
%==========================================================================

Eval=diag(Eval);
[~,ind]=sort(Eval,'descend'); % Index of growing eigenvalue
Eval=Eval(ind);
Evec=Evec(:,ind);


u1=Evec(:,1);
u1=u1./norm(u1); % Normalised eigenvector in growing direction (used for stable manifold)

u2=Evec(:,2);
u2=u2./norm(u2); % Normalised eigenvector in decaying direction (used for unstable manifold)

%==========================================================================
% Stable manifold
%==========================================================================

q_0=1;
p_0=1;

% One side

q_stable=q_fp-q_0*mu*u1(1);
p_stable=p_fp-p_0*mu*u1(2);
z=[q_stable,p_stable,0];
[t,z_ps] = ode89(@(t,z) dhat(t,z,omega,beta,delta,gamma,system_flag,time_flag) ,[0 t_final],z); % Integrate
    
figure(1)
hold on
plot(z_ps(:,1),z_ps(:,2),'r-','Markersize',1)
xlabel('q')
ylabel('p')
axis([qmin qmax pmin pmax])

% Other side

q_stable=q_fp+q_0*mu*u1(1);
p_stable=p_fp+p_0*mu*u1(2);
z=[q_stable,p_stable,0];
[~,z_ps] = ode89(@(t,z) dhat(t,z,omega,beta,delta,gamma,system_flag,time_flag) ,[0 t_final],z); % Integrate
   
figure(1)
hold on
l2=plot(z_ps(:,1),z_ps(:,2),'r-','Markersize',1);
xlabel('q')
ylabel('p')
axis([qmin qmax pmin pmax])

%==========================================================================
% Unstable manifold
%==========================================================================

time_flag='bwd'; % Time-reversed dynamics
t_final=t_final_bwd;

% One side
q_unstable=q_fp+q_0*mu*u2(1); % Now we use the 
p_unstable=p_fp+p_0*mu*u2(2);
z=[q_unstable,p_unstable,0];
[~,z_ps] = ode89(@(t,z) dhat(t,z,omega,beta,delta,gamma,system_flag,time_flag) ,[0 t_final],z); % Integrate
        
 figure(1)
hold on
plot(z_ps(:,1),z_ps(:,2),'b-','Markersize',1)
xlabel('q')
ylabel('p')
axis([qmin qmax pmin pmax])        

% Otherside

q_unstable=q_fp-q_0*mu*u2(1);
p_unstable=p_fp-p_0*mu*u2(2);
z=[q_unstable,p_unstable,0];
[~,z_ps] = ode89(@(t,z) dhat(t,z,omega,beta,delta,gamma,system_flag,time_flag) ,[0 t_final],z); % Integrate
        
        
figure(1)
hold on
l3=plot(z_ps(:,1),z_ps(:,2),'b-','Markersize',1);
xlabel('q')
ylabel('p')
legend([l1 l2 l3],{'Fixed point','Stable manifold','Unstable manifold'})
axis([qmin qmax pmin pmax])
box on
pause(1)
       
    




