function [J_out,psi,E]=get_jacobian_mexican_hat(q,p,omega,beta,delta,gamma,system_flag) 


bad_system_flag=0; % Counters to check if bad system flag passed

%==========================================================================
% Jacobian of the Hamitlonian K 
%==========================================================================

if ismember(system_flag,'K')

bad_system_flag=1; 


J_qq=2*beta*q*p-gamma;
J_pq=-omega+beta*q^2+3*beta*p^2;
J_qp=omega-beta*p^2-3*beta*q^2;
J_pp=-2*beta*q*p-gamma;
J_out=[J_qq J_pq; J_qp J_pp];
[psi, E] = eigs(J_out);

end

if ismember(system_flag,'MKS')

bad_system_flag=1; 


J_qq=-2*beta*q*p-gamma;
J_pq=omega-beta*q^2-3*beta*p^2;
J_qp=-omega+beta*p^2+3*beta*q^2;
J_pp=2*beta*q*p-gamma;
J_out=[J_qq J_pq; J_qp J_pp];
[psi, E] = eigs(J_out);

end

%==========================================================================
% Checks
%==========================================================================

if bad_system_flag==0  
    'Bad system flag'
    system_flag
    return
end

end



