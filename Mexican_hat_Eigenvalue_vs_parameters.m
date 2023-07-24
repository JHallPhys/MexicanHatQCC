%Uses code Mexican_hat_Eigenvalue_number
%Uses parfor - parallel computing package: 

tic 
D=100; %200 takes 3 mins. 300 takes 5 mins. 400 is 9mins

gamma = linspace(0.00,1,D);
delta = linspace(0.00,10,D);

Num_evalues_a = zeros(length(gamma),length(delta));
Num_evalues_b = zeros(length(gamma),length(delta));
Num_evalues_c = zeros(length(gamma),length(delta));
Num_branches = zeros(length(gamma),length(delta));

for j = 1:length(gamma)
    parfor n = 1:length(delta)
[Number_evalues_a,Number_evalues_b,Number_evalues_c,branches] = Mexican_hat_Eigenvalue_number(gamma(j), delta(n));
 %a is the "straight line"/degenerate. b is the curve connected but a that
 %goes upwards (contains the ground state) and c is the infinite curve that
 %goes downward. 
    Num_evalues_a(n,j) = Number_evalues_a;
    Num_evalues_b(n,j) = Number_evalues_b;
    Num_evalues_c(n,j) = Number_evalues_c;
    Num_branches(n,j) = branches;
    end
%     if j == mod(10)
%         display('1/10 done')
%     end
end
whos

figure(1);
surf(gamma,delta,Num_evalues_a,'linestyle','none')
title('Number of Eigenvalues in left-most section')
xlabel('Gamma') 
ylabel('Delta')
%caxis([0 max(number)])
view(2)

figure(2);
surf(gamma,delta,Num_evalues_b,'linestyle','none')
title('Number of Eigenvalues in top curve')
xlabel('Gamma') 
ylabel('Delta') 
view(2)

figure(3);
surf(gamma,delta,Num_evalues_c,'linestyle','none')
title('Number of Eigenvalues in bottom curve')
xlabel('Gamma') 
ylabel('Delta') 
view(2)

figure(4);
surf(gamma,delta,Num_branches,'linestyle','none')
title('Number of eigenvalue branches')
xlabel('Gamma') 
ylabel('Delta') 
view(2)

toc
