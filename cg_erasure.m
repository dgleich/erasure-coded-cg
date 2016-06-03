function [x, flag, iter, fail_ind, x_fail_ind, resvec] = cg_erasure(A, b, tol, maxit, x0, fail_ind, verbose, varargin)
% FUNCTION cg_erasure is conjugate gradient for solving a sqaure linear system
%                       Ax = b          (1)
% which is used to experiment with the idea of erasure coded computation.
% For a linear system defined in (1), the augmented system through erasure
% coding the rows and columns of matrix A is defined as
%               [A A*E; E'*A E'*A*E][x; zeros] = [b; E'*b] (2)
% where E is the encoding matrix.
% Thus from the augmented system point of view, we chase at the unique
% solution to (2) with the last component being 0. This is feasible because
% we know the null basis of the augmented matrix is [e;-1]. To recover the
% solution x, we can do the following:
%   - for all non-fail components, just add the value of the augmented
%   component to it;
%   - for the only fail component, either we retrieve it by substituting
%   the rest (n-1) to one equation of (1), or we use x_fail_ind, which is
%   its value just before halting, and add to it the value of augmented
%   component.
% A, the coefficient matrix, or a function handle, possibly after erasure
% coding.
% b, the RHS, possibly after erasure coding.
% tol, stopping criterion in relative residual.
% maxit, maximum number of iterations.
% x0, initial guess.
% fail_ind, the index of the failure component. If it is non-empty, we
% assume the failure happens before the CG starts.
% verbose, a flag indicates whether print out the relative residual of each
% iteration.
% varargin{1}, the value of the failure point.
% return - x, approximate solution.
% return - flag, indication of whether convergence is met.
%          flag = 0, converged,
%          flag = 1, not converged.
% return - iter, the actual number of iterations used.
% return - fail_ind, the index of the failure component.
% return - x_fail_ind, the failed component of x just before halting.
% return - resvec, the vector of residual history.
% History: 2014-08-27. I add an input parameter fail_ind, which is the
%          index of the failed component. We require 1<=fail_ind<=(n-1),
%          where the n-th row and column of A is for erasure code. To
%          simulate failure, the code is written in the way that before and
%          after any aggregation operation, the failed component is zeroed
%          out.
%          2014-08-28. I add the input parameter varargin, where
%          varargin{1} is the iteration number at which the failure
%          happens. For simplicity, we assume the failure can only happen
%          after the varargin{1}-th iteration completes. Admittedly, a
%          meaningful value of varargin{1} can only be determined after we
%          know how many iterations CG need when running on the augmented
%          system when NO failure happens.
%          2014-09-07. I extend the code to allow for multiple failure
%          components. For this purpose, I let varargin{2} be the number of
%          failure components. If varargin{2} is not given, we let the
%          default number of failure components be 1.
%          2014-10-12. I modified the code in such a way that no aggregation
%          operation (espcially the inner-product) result is reused. We
%          also adde the code to record the residual norm history.

if nargin > 7
    fail_point = varargin{1};
else
    fail_point = [];
end

if nargin > 8
    num_fail = varargin{2};
else
    num_fail = 1;
end

flag = 1; % assume no convergence.

n = length(b);

if isempty(x0)
    x = zeros(n,1);
else
    x = x0;
end

% fix b_nrm at 1 if we want to use absolute residual.
%b_nrm = 1;

% before compute b_nrm, we zero out the failed component.
if ~isempty(fail_ind)
    b(fail_ind) = 0;
end
b_nrm = norm(b);
if b_nrm < 100*eps
    b_nrm = 1;
end

% before act on x, we zero out the failed component.
if ~isempty(fail_ind)
    x(fail_ind) = 0;
end
if isa(A, 'function_handle')
    Ax = A(x);
else
    Ax = A*x;
end
% after matvec, we zero out the failed component.
if ~isempty(fail_ind)
    Ax(fail_ind) = 0;
end

r = b - Ax;
% before compute norm(r), we zero out the failed component.
if ~isempty(fail_ind)
    r(fail_ind) = 0;
end
%if norm(r)/b_nrm < tol
if norm(r) < tol
    flag = 0;
    iter = 0;
    return;
end

x_fail_ind = [];

% this flag is used to inform the update to the descent direction p whether
% any fault happen in between two iterations of CG.
fault_flag = 0;

resvec = zeros(maxit+1,1);

for iter = 1:maxit
 % before compute the inner product, we zero out the failed component.
 if ~isempty(fail_ind)
     r(fail_ind) = 0;
 end
 rho = (r'*r);
 % record the norm of the initial residual.
 if iter == 1
     resvec(1) = norm(r);
 end
 
 if iter > 1
    % Yao Zhu added on 2014-10-12. Recompute rho_old simultaneously with rho.
    % before compute the inner product, we zero out the failed component.
    if ~isempty(fail_ind)
     r_old(fail_ind) = 0;
    end
    rho_old = (r_old'*r_old);
    
    beta = rho / rho_old;
    % check whether any fault has happened.
    if ~fault_flag
        p = r + beta*p;
    else
        % Yao Zhu added on 2014-10-13.  maintain r_old in p still does not
        % work. I think the main reason is that r and r_old are not
        % mutually orthogonal when restricted on the degenerated subsystem.
        % for this purpose, we should just use the restart strategy. That
        % is making p=r.
        %p = r + beta*r_old;     % truncate.
        p = r;     % truncate.
        fault_flag = 0;         % indicate truncation of p already done.
    end
 else
    p = r;
 end

 % before act on p, we zero out the failed component.
 if ~isempty(fail_ind)
     p(fail_ind) = 0;
 end
 % This MatVec could have one component failed.
 if isa(A, 'function_handle')
     q = A(p);
 else
     q = A*p;
 end
 % after matvec, we zero out the failed component.
 if ~isempty(fail_ind)
     q(fail_ind) = 0;
 end

 % Yao Zhu added on 2014-10-12. Recompute rho simultaneously with (p'*q).
 % before compute the inner product, we zero out the failed component.
 if ~isempty(fail_ind)
  r(fail_ind) = 0;
 end
 rho = (r'*r);
 
 % before compute the inner product, we zero out the failed component.
 if ~isempty(fail_ind)
     p(fail_ind) = 0;
     q(fail_ind) = 0;
 end
 alpha = rho / (p'*q);
 % update approximate solution.
 x = x + alpha * p;
 
 % Yao Zhu added on 2014-10-12. record r_old for the purose of recomputing
 % rho_old.
 r_old = r;
 
 % update residual.
 r = r - alpha*q;
 
 % check convergence 
 % before compute norm(b), we zero out the failed component.
 if ~isempty(fail_ind)
    r(fail_ind) = 0;
    % when failure happens, we also need to recompute b_nrm after zero out
    % the failed component.
    b(fail_ind) = 0;
    b_nrm = norm(b);
    if b_nrm < 100*eps
        b_nrm = 1;
    end
 end
 
 if verbose
     %fprintf('cg_erasure iteration: %d, relative residual: %e\n', iter, norm(r) / b_nrm); 
     fprintf('cg_erasure iteration: %d, norm of residual: %e\n', iter, norm(r));
     resvec(iter+1) = norm(r);
 end
 
 % use relative residual.
 %if norm(r) / b_nrm <= tol
 % use norm of residual.
 if norm(r) <= tol
     flag = 0;
     break;
 end
 
 % For simplicity, we only assume failure can happen after an iteration
 % completes. And this intermediate failure can happen only when fail_ind
 % is empty. Then we introduce a failed component in the first (n-1)
 % components.
 if ~isempty(fail_point) && isempty(fail_ind) && (num_fail > 0)
    if iter == fail_point
        fail_ind = randperm(n-num_fail);
        fail_ind = fail_ind(1:num_fail);
        % check the value of x(fail_ind) just before it halts.
        %disp(['just before halting, x(fail_ind)= ']);
        %disp(x(fail_ind));
        x_fail_ind = x(fail_ind);
        % Yao Zhu added the following on 2014-09-26. The motivation is that
        % we forget the information of failure components totally.
        % now we need to forget about components x(fail_ind).
        %x(fail_ind) = 0;
        %%%%x=zeros(size(x));
        % We also need to delete the information of x(fail_ind) from r. we
        % also need to forget the conjugate direction p totally.
        %r = b - A*x;
        
        % Yao Zhu added on 2014-10-13. When failure happens, it's very
        % important to truncate the expansion of p in r (for the expression
        % see page 53 of (Meurant 2006).
        fault_flag = 1;   % record fault happened before next iteration.
        %%%%p = zeros(size(p));  % we now use fault_flag to inform the
                                 % statement update p, instead of
                                 % truncating here.
    end
 end
 
end % for iter = 1:maxit

resvec(iter+2:end) = [];