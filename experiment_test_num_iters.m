%% The goal of this experiment is to test the number of iterations of CG
% with and without the erasure coding to see how it changes if we add some
% number of columns 

% A lot of this comes from run_script_cg_erasure.m

matrix_file = 'Ltridiag500';

disp(['matrix_file is: ' matrix_file]);
load([matrix_file '.mat']);
n = size(A,1);

x_true=rand(n,1);
b =  A*x_true;


tol = 1e-10;
n = size(A, 1);
maxit = 10*n; % the maximum number of cg iterations.

%%

num_fail = 0;

E = zeros(n,num_fail);
x0 = [];
fail_ind = [];          % start with no failure.
fail_point = 0;

A2=[A A*E; E'*A, E'*A*E]; % the erasure coded augment matrix.
x2=[x_true; zeros(num_fail,1)];     % fix the redundant elements at 0 in the augmented solution.
b2=A2*x2;           % note b2 is nothing but [b; E'*b]

[x2_sol_fail, flag, iter, fail_ind, x_fail_ind, resvec] = cg_erasure(A2, b2, tol, maxit, x0, fail_ind, 1, fail_point, num_fail);

%%
num_fail = 10;

% we sample the encoding matrix from Gaussian ensemble.
%randn('seed', 0);
E = randn(n, num_fail)/(sqrt(n)); % properly scaling E is very important to the condition number of A2.

A2=[A A*E; E'*A, E'*A*E]; % the erasure coded augment matrix.
x2=[x_true; zeros(num_fail,1)];     % fix the redundant elements at 0 in the augmented solution.
b2=A2*x2;           % note b2 is nothing but [b; E'*b]
% the null basis of A2:
null_basis = [E; -eye(num_fail)];

% Report condition number
lam1 = eig(full(A));
lam2 = eig(full(A2)); 
lam2(lam2 < n*eps(max(lam2))) = []; % remove null-space
disp(['Condition number of  A = ', num2str(max(lam1)/min(lam1))]);
disp(['Condition number of At = ', num2str(max(lam2)/min(lam2))]);

x0 = [];
fail_ind = [];          % start with no failure.
fail_point = 0;
[x2_sol_fail, flag, iter, fail_ind, x_fail_ind, resvec] = cg_erasure(A2, b2, tol, maxit, x0, fail_ind, 1, fail_point, num_fail);

%disp(['Number of cg_erasure iterations = ' num2str(iter)]);

%%
