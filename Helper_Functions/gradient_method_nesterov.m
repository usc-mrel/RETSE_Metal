function [Tx,costOut] = gradient_method_nesterov(I1,I2,Tx,Btotal,Bgrad,...
                                                    x0,opts)
% INPUT
% ======================================
% I1, I2..  Distorted images
% Tx .....  Current displacement estimate
% Btotal..  Sparse GPU array containing the transform matrix from 
%           spline coefficients to spatial distortion.
% Bgrad...  Sparse GPU array containing the gradient of the transform
%           matrix in distorted dimension
% x0......  Initial guess of B-spline coefficients
% opts ...  Contains stepsize and termination criteria

% ======================================
% OUTPUT
% ======================================
% Tx .....  The optimal solution within specified tolerance
% cost ..   The value of the objective function at xopt

maxIter = opts.maxIter;
epsilon = opts.epsilon;
stepsize = opts.stepsize;
stag_lim = opts.stag_lim;
stag_tol = opts.stag_tol;

x = x0;
y = x0;
iter = 0;
lambda = 0;

grad = gradstep(I1,I2,Tx,Btotal,Bgrad);

[nx, ny, nz] = size(Tx);

costOut = zeros(maxIter,1);

stag_ctr = 0;
costOld = Inf;

while norm(grad(:)) > epsilon && iter < maxIter && ...
        stag_ctr < stag_lim
    iter = iter + 1;
    
    lambdaNew = (1 + sqrt(1 + 4*lambda^2))/2;
    gamma = (1 - lambda)/lambdaNew;
    
    yNew = x - stepsize * grad;
    xNew = (1 - gamma) * yNew + gamma * y;
    
    x = xNew;
    y = yNew;
    lambda = lambdaNew;
    
    Tx = reshape((xNew(:).' * Btotal).',[nx, ny, nz]);
    grad = gradstep(I1,I2,Tx,Btotal,Bgrad);
    cost = cost_func(I1,I2,Tx);
    
    costOut(iter) = gather(cost);
    
    if abs(costOut(iter) - costOld)/costOld < stag_tol
        stag_ctr = stag_ctr + 1;
    else 
        stag_ctr = 0;
    end
    
    costOld = costOut(iter);
    
    fprintf('iter = %d \t norm(grad) = %2.6f \t cost = %2.6f \n',...
        iter,norm(grad),cost);
end

