function [Number_evalues_a,Number_evalues_b,Number_evalues_c,branches] = Mexican_hat_Eigenvalue_number(gamma, delta)
N=100; %matrix size

%matrices for postition and momentum operators
n=1:N-1;
k=sqrt(n);

a = diag(k,1);
adagger = diag(k,-1);

%diag(k,1) is a and diag(k,-1) is adagger
q=sqrt(0.5)*(adagger + a);
p=1i*sqrt(0.5)*(adagger - a);

%Mexican hat
omega = -1.0;
gamma = gamma;
beta = 0.05;
zeta = 0.0;
delta = delta;
epsilon=0;
%H = 0.5*(omega-1i*gamma)*(p^2+q^2) +0.25*(beta-1i*delta)*(p^2+q^2)^2 ...
%    + (kappa-1i*epsilon)*q;

H = (omega-1i*gamma)*(adagger*a) +(beta-1i*zeta)*adagger*adagger*a*a ...
    + (delta-1i*epsilon)*(a+adagger)/sqrt(2); %+ 1/2*eye(N);

[eigenvectors, eigenvalues] = eig(H);

[V,E] = sortem( eigenvectors,eigenvalues);

eigenvalues_im = E; %Ordered by imaginary part

eigenvalues_r = sort(eigenvalues_im,'ComparisonMethod','real'); %ordered by real part

eigenvalues_a = zeros(1,length(eigenvalues_r));
eigenvalues_b = zeros(1,length(eigenvalues_r));
eigenvalues_c = zeros(1,length(eigenvalues_r));

%   if delta ~=0
% eigenvalues_a(1) = eigenvalues_r(1);
% for n = 2:length(eigenvalues_r)
%   eigenvalues_a(n) = eigenvalues_r(n); %This line might need tweaking for when there is only one or no values in a
%   if imag(eigenvalues_r(n))<imag(eigenvalues_r(n-1))
%       for m = n:length(eigenvalues_r)
%           if imag(eigenvalues_r(m))<imag(eigenvalues_r(n-1))
%               eigenvalues_c(m) = eigenvalues_r(m);
%           else
%               eigenvalues_b(m) = eigenvalues_r(m);
%           end
%       end
%       break
%   end
% end
%   end

if delta ~=0
eigenvalues_a(1) = eigenvalues_r(1);
for n = 2:length(eigenvalues_r)
  if imag(eigenvalues_r(n))<imag(eigenvalues_r(n-1))
      parfor m = n:length(eigenvalues_r)
          if imag(eigenvalues_r(m))<=imag(eigenvalues_r(n-1))
              eigenvalues_c(m) = eigenvalues_r(m);
          else
              eigenvalues_b(m) = eigenvalues_r(m);
          end
      end
      break
  end
  eigenvalues_a(n) = eigenvalues_r(n);
end
end


if delta == 0
      eigenvalues_b(1) = eigenvalues_im(1);
    parfor n = 2:length(eigenvalues_im)
        if real(eigenvalues_im(n-1))<real(eigenvalues_im(n))
        eigenvalues_b(n) = eigenvalues_im(n);
        else
        eigenvalues_c(n) = eigenvalues_im(n);
        end
    end
end 
%eigenvalues_c %7.51 - 7.12i
%eigenvalues_b(26)   %8.7841 - 7.0879i
%9.68 - 7.36i



Number_evalues_a = nnz(eigenvalues_a);
Number_evalues_b = nnz(eigenvalues_b);
Number_evalues_c = nnz(eigenvalues_c);

branches = 3;
if any(Number_evalues_a) == 0
branches = branches -1;
end
if any(Number_evalues_b) == 0
branches = branches -1;
end
if any(Number_evalues_c) == 0
branches = branches -1;
end

end
