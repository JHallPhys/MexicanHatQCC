function [qfp_out,pfp_out]=get_fixed_points_mexican_hat(omega,beta,delta,gamma,system_flag) 

% Want to take in the parametrs as real harmonic, real quartic, tilt etc
bad_system_flag=0; % Counters to check if bad system flag passed

%==========================================================================
% Fixed Points of the Hamiltinian K
%==========================================================================

if ismember(system_flag,'K')
bad_system_flag=1;
% Define the symbolic variables
syms qfp pfp

% Define the equations
qdot=-omega*pfp+beta*(qfp^2+pfp^2)*pfp-gamma*qfp;
pdot=omega*qfp-beta*(pfp^2+qfp^2)*qfp-delta-gamma*pfp;
% Solve the system of equations
solutions = solve([qdot, pdot], [qfp, pfp]);

% Display the roots
disp(" ")
disp("-----------");
disp("Fixed Points of the system:");
disp("-----------");
for i = 1:numel(solutions.qfp)
    fprintf('x%d = %f\n', i, solutions.qfp(i));
    fprintf('y%d = %f\n', i, solutions.pfp(i));
    disp("-----------");
end

qfp_out=double(vpa(solutions.qfp));
pfp_out=double(vpa(solutions.pfp));


end

%==========================================================================
% Fixed Points of the Hamiltinian -K*
%==========================================================================

if ismember(system_flag,'MKS')
bad_system_flag=1;
% Define the symbolic variables
syms qfp pfp

% Define the equations
qdot=omega*pfp-beta*(qfp^2+pfp^2)*pfp-gamma*qfp;
pdot=-omega*qfp+beta*(pfp^2+qfp^2)*qfp+delta-gamma*pfp;
% Solve the system of equations
solutions = solve([qdot, pdot], [qfp, pfp]);

% Display the roots
disp(" ")
disp("-----------");
disp("Fixed Points of the system:");
disp("-----------");
for i = 1:numel(solutions.qfp)
    fprintf('x%d = %f\n', i, solutions.qfp(i));
    fprintf('y%d = %f\n', i, solutions.pfp(i));
    disp("-----------");
end

qfp_out=double(vpa(solutions.qfp));
pfp_out=double(vpa(solutions.pfp));


end


%==========================================================================
%Checkz
%==========================================================================


if bad_system_flag==0  
    'Bad system flag'
    system_flag
    return
end

end