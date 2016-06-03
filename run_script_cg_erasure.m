% script to run erasure coded cg.
% we need to set the following variables manually:
% matrix_file - the file storing the SPD raw matrix A.
% tol - the stopping criterion for CG.
% num_fail - the number of failure components.

disp(['matrix_file is: ' matrix_file]);
load([matrix_file '.mat']);

n = size(A, 1);
maxit = 10*n; % the maximum number of cg iterations.
% a fractional num_fail is interpreted as fault rate.
if num_fail < 1 && num_fail > 0
    num_fail = ceil(n*num_fail);
end

% for sanity check, we'd better use a random x.
%rand('seed', 0);
x_true=rand(n,1);
% make up the rhs.
b =  A*x_true;

% we sample the encoding matrix from Gaussian ensemble.
%randn('seed', 0);
E = randn(n, num_fail)/sqrt(n); % properly scaling E is very important to the condition number of A2.

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
% sample the failure point from within 0.25*n.
if num_fail
    fail_point = randperm(ceil(0.25*n),1);
else 
    fail_point = 0;
end



% x2_sol_fail, the solution with cg_erasure continuing after a component failed.
% flag, indicates whether cg_erasure converged.
% iter, the number of iterations used.
% fail_ind, index of the failed component/processor (assume each processor
% for a row).
% x_fail_ind, the logged value of the failed component of x just before the
% failure happens.
[x2_sol_fail, flag, iter, fail_ind, x_fail_ind, resvec] = cg_erasure(A2, b2, tol, maxit, x0, fail_ind, 1, fail_point, num_fail);

disp(['Number of cg_erasure iterations = ' num2str(iter)]);

%when cg_erasure terminates, we first check whether for the subsystem
%defined on the non-failure components, it gives the solution satisfying
%the tol.
if ~isempty(fail_ind)
    A2_sub = A2;
    A2_sub(fail_ind, :) = [];
    A2_sub(:, fail_ind) = [];
    b2_sub = b2;
    b2_sub(fail_ind) = [];
    x2_sol_fail_sub = x2_sol_fail;
    x2_sol_fail_sub(fail_ind) = [];
    %disp('relative residual on the non-failure component subsystem');
    %norm(b2_sub-A2_sub*x2_sol_fail_sub) / norm(b2_sub)
    % Yao Zhu added on 2014-10-13. Note the above relative residual on the
    % non-failure component subsytem could be large. This is because b2_sub
    % is NOT the correct purified RHS from incoporating the snapshot values
    % of the failed components.
end

%when cg_erasure terminates, it's very important to check the quality on
%the global system (A2, b2) as achieved by x2_sol_fail.
%disp(['relative residual on the encode linear system is: ' num2str(norm(b2-A2*x2_sol_fail) / norm(b2))]);
disp(['Residual norm on the encode linear system is: ' num2str(norm(b2-A2*x2_sol_fail))]);
disp(['Compare with the tol for cg_erasure we intend: ' num2str(tol)]);

% check whether x2_sol_fail(fail_ind) is fixed at x_fail_ind
if ~isempty(fail_ind)
    disp('Check whether x2_sol_fail(fail_ind) is fixed at x_fail_ind...');
    disp(['norm(x2_sol_fail(fail_ind)-x_fail_ind) = ' num2str(norm(x2_sol_fail(fail_ind)-x_fail_ind))]);
end

% recover the intended solution.
disp('Recover the intended solution...');
x2_sol_fail = x2_sol_fail + null_basis * x2_sol_fail(end-num_fail+1:end);

% check the difference with the intended solution.
disp('Check the difference with the intended solution...');
disp(['On the encode solution norm(x2 - x2_sol_fail)/norm(x2) = ' num2str(norm(x2 - x2_sol_fail)/norm(x2))]);
disp(['On the raw solution norm(x_true - x2_sol_fail(1:n))/norm(x_true) = ' num2str(norm(x_true - x2_sol_fail(1:n))/norm(x_true))]);

% it makes more sense to check the relative residual on the raw system than
% check the relative error wrt to the intended solution.
disp('Check relative residual on the raw system...');
disp(['On the raw solution norm(b-A*x2_sol_fail(1:n))/norm(b) = ' num2str(norm(b-A*x2_sol_fail(1:n))/norm(b))]);