function [Tx,cost] = calc_bfield(I1,I2,varargin)

% INPUT
% ======================================
% I1, I2..  Distorted images, with the first dimension being distorted
% varargin  Contains stepsize and termination criteria

% ======================================
% OUTPUT
% ======================================
% Tx .....  The optimal solution within specified tolerance
% cost ..   The value of the objective function at xopt
    nin = abs(nargin) - 2;

    if nin >= 1
        opts = varargin{1};
    else
        opts = struct();
    end

    if ~isstruct(opts)
        fprintf('Optional input not a struct. Resorting to defaults. \n');
    else

        if ~isfield(opts,'maxIter')
            opts.maxIter = 400;
        end

        if ~isfield(opts,'stag_tol')
            opts.stag_tol = 6.5e-4;
        end

        if ~isfield(opts,'stag_lim')
            opts.stag_lim = 10;
        end

        if ~isfield(opts,'stepsize')
            opts.stepsize = 0.1;
        end

        if ~isfield(opts,'epsilon')
            opts.epsilon = 1e-4;
        end
        
        if ~isfield(opts,'spacing')
            opts.spacing =[2,2,2];
        end
    end

    I1 = gpuArray(I1);
    I2 = gpuArray(I2);

    [nx, ny, nz] = size(I1);

    Tx = zeros(nx,ny,nz,'gpuArray');

    %these are in voxel units.
    hx = opts.spacing(1);
    hy = opts.spacing(2);
    hz = opts.spacing(3);

    %only ever need the transposes, so work with them directly.
    Bx = sparse(bsplines(nx,hx)).';
    By = sparse(bsplines(ny,hy)).';
    Bz = sparse(bsplines(nz,hz)).';

    %number of splines in each dimension.
    mx = size(Bx,1);
    my = size(By,1);
    mz = size(Bz,1);

    Btotal = gpuArray(sparse(kron(Bz,kron(By,Bx))));

    Bgx = sparse(gradient(Bx));
    Bgrad = gpuArray(sparse(kron(Bz,kron(By,Bgx))));

    x0 = zeros(mx*my*mz,1);

    [Tx,cost] = gradient_method_nesterov(I1,I2,Tx,Btotal,Bgrad,x0,opts);

    fprintf('Finished \n');
end




