% This script tests the conditioning of the purified linear system. This
% corresponds to the system actually solved by the fault-tolerance
% encoding. N
matrix_file = 'Ltridiag500';
matrix_file = 'mhdb416';
matrix_file = 'nos3';

disp(['matrix_file is: ' matrix_file]);
load([matrix_file '.mat']);

n = size(A,1);

ntrials = 10;
num_fails = [1,2,5,10,25,50,100];

for fi=1:numel(num_fails)
    num_fail = num_fails(fi);
    disp(' ');
    disp(['Number of failures = ', num2str(num_fail)]);
    disp(' ');
    for ti=1:ntrials

        % we sample the encoding matrix from Gaussian ensemble.
        %randn('seed', 0);
        E = randn(n, num_fail)/sqrt(n); % properly scaling E is very important to the condition number of A2.

        A2=[A A*E; E'*A, E'*A*E]; % the erasure coded augment matrix.
        %x2=[x_true; zeros(num_fail,1)];     % fix the redundant elements at 0 in the augmented solution.
        % the null basis of A2:
        null_basis = [E; -eye(num_fail)];

        % Report condition number
        symmat = @(A) triu(A,1) + triu(A,1)' + diag(diag(A));
        lam1 = eig(full(A));
        lam2 = eig(full(symmat(A2))); 
        lam2(lam2 < n*eps(max(lam2))) = []; % remove null-space
        disp(['Condition number of  A = ', num2str(max(lam1)/min(lam1))]);
        disp(['Condition number of At = ', num2str(max(lam2)/min(lam2))]);

    end    
end    