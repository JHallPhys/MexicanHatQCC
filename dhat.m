function pdr=dhat(t,y,omega,beta,delta,gamma,system_flag,time_flag) 

bad_time_flag=0; % Counters to check if bad time flag passed
bad_system_flag=0; % Counters to check if bad system flag passed


%==========================================================================
% Equations of motion for mexican hat Hamitlonian
%==========================================================================

if ismember(system_flag,'K')

    bad_system_flag=1;
    q=y(1);
    p=y(2);
    
    if ismember(time_flag,'fwd') % Forward in Time 
         bad_time_flag=1;
    elseif ismember(time_flag,'bwd') % Backward in time
         bad_time_flag=1; 
         
         omega=-omega;
         beta=-beta;
         delta=-delta;
         gamma=-gamma;
    end

    qdot=-omega*p+beta*(q^2+p^2)*p-gamma*q;
    pdot=omega*q-beta*(p^2+q^2)*q-delta-gamma*p;
    pdr=[qdot;pdot;-2*gamma*(q^2+p^2)*y(3)];

end

%==========================================================================
% Equations of motion for mexican hat Hamitlonian
%==========================================================================

if ismember(system_flag,'MKS')

    bad_system_flag=1;
    q=y(1);
    p=y(2);
    
    if ismember(time_flag,'fwd') % Forward in Time 
         bad_time_flag=1;
    elseif ismember(time_flag,'bwd') % Backward in time
         bad_time_flag=1; 
         
         omega=-omega;
         beta=-beta;
         delta=-delta;
         gamma=-gamma;
    end

    qdot=omega*p-beta*(q^2+p^2)*p-gamma*q;
    pdot=-omega*q+beta*(p^2+q^2)*q+delta-gamma*p;
    pdr=[qdot;pdot;-2*gamma*(q^2+p^2)*y(3)];

end



%==========================================================================
% Idiot proofing: fail on value 0 
%==========================================================================



if bad_time_flag==0  
    'Bad time flag'
    time_flag
    return
end

if bad_system_flag==0  
    'Bad system flag'
    system_flag
    return
end


end