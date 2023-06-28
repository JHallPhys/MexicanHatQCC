clear all 
close all




%==========================================================================
% Initialisation
%==========================================================================

N=11; % Number of initial conditions

omega=1;%
gamma=0.05
beta=0.05;
delta=1;
t_final=20;

qmin=-10
qmax=10
pmin=-10
pmax=10

% qmin=-2
% qmax=2
% pmin=-2
% pmax=2

q=linspace(qmin,qmax,N );
p=linspace(pmin,pmax,N);
system_flag='MKS';
time_flag='fwd';
%==========================================================================
% Get Trajectories 
%==========================================================================


reverseStr = '';

for j=1:length(q)

    % Counter for the number we are on
    msg = sprintf('GetNorm %d/%d', j, length(q));
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));


    for k = 1:length(p)
        
        % 1) Initial conditions in phase space       
        z(1)=q(j);
        z(2)=p(k);
        z(3)=1;
        
          % 2) Propagate Norm
        
        [t,z_ps] = ode89(@(t,z) dhat(t,z,omega,beta,delta,gamma,system_flag,time_flag) ,[0 t_final],z); % Integrate
        
        
    figure(1)
    hold on
    plot(z_ps(:,1),z_ps(:,2),'k-','Markersize',1)
    xlabel('q')
    ylabel('p')
    axis([qmin qmax pmin pmax])
    box on


        
    end
       
    
end
% return
%==========================================================================
% Plot the fixed points
%==========================================================================

[q_fp,p_fp]=get_fixed_points_mexican_hat(omega,beta,delta,gamma,system_flag); 

figure(1)
hold on 
plot(q_fp(1),p_fp(1),'r.','Markersize',10)
plot(q_fp(2),p_fp(2),'b.','Markersize',10)
plot(q_fp(3),p_fp(3),'r.','Markersize',10)
