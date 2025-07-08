function outp_dual = type1_pec_sense_normal_dual(inp, sens, mask, nj, p1, p2, iflag, eps, support_mask, lambda, topup_sign, phase)

%% Start a stopwatch timer
start_time = tic;

%% Declare a persistent variable
persistent cg_iter;
if isempty(cg_iter)
    cg_iter = 0;
else
    cg_iter = cg_iter + 1;
end

%% Get imaging parameters
[Nkx,Nky,Nkz,Nc,Nacq] = size(sens);

%--------------------------------------------------------------------------
% Calculate the total number of voxels
%--------------------------------------------------------------------------
N = Nkx * Nky * Nkz;

outp_dual = zeros(N,1);

for nacq = 1:Nacq
    %--------------------------------------------------------------------------
    % Calculate the number of voxels within a support mask
    %--------------------------------------------------------------------------
    N_support = length(find(support_mask(:,:,nacq)));

    %% Calculate scale factors
    type1_scale_factor = 1 / sqrt(Nkx * Nky);
    type2_scale_factor = 1 / sqrt(Nkx * Nky);

    %% Calculate S * m
    inp_md = reshape(inp, [Nkx Nky Nkz]);
    Sm = bsxfun(@times, sens(:,:,:,:,nacq), exp(1i*phase(:,:,nacq)) .* inp_md); % Nkx x Nky x Nkz x Nc
    Sm = reshape(Sm, [N Nc]);

    %% Calculate (I_{Nc} kron F) * S * m
    %--------------------------------------------------------------------------
    % Perform type-1 NUFFT (nonuniform, image space => uniform, k-space)
    % Siemens: k-space <=> image space
    %--------------------------------------------------------------------------
    tstart = tic; fprintf('%s:(PCG=%d) Performing type-1 NUFFT... ', datetime, cg_iter);
    cj = reshape(double(Sm((support_mask(:,:,nacq) > 0),:)), [N_support Nc]);
    FSm = reshape(single(type1_scale_factor * finufft2d1(p1{nacq}, p2{nacq}, cj, -iflag, eps, Nkx, Nky)), [Nkx Nky Nkz Nc]);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %% Calculate d = (I_{Nc} kron R_{Omega}) * (I_{Nc} kron F) * S * m
    d = bsxfun(@times, FSm, mask(:,:,nacq)); % Nkx x Nky x Nkz x Nc

    %% Calculate (I_{Nc} kron R_{Omega}^H) * d
    Rhd = bsxfun(@times, mask(:,:,nacq), d);

    %% Calculate (I_{Nc} kron F^H) * (I_{Nc} kron R_{Omega}^H) * d
    %--------------------------------------------------------------------------
    % Perform type-2 NUFFT (nonuniform, image space <= uniform, k-space)
    % Siemens: k-space <=> image space
    %--------------------------------------------------------------------------
    tstart = tic; fprintf('%s:(PCG=%d) Performing type-2 NUFFT... ', datetime, cg_iter);
    FhRhd = complex(zeros(N, Nc, 'single'));
    cj = type2_scale_factor * finufft2d2(p1{nacq}, p2{nacq}, iflag, eps, reshape(double(Rhd(:,:,1,:)), [Nkx Nky Nc])); % N x Nc
    FhRhd((support_mask(:,:,nacq) > 0),:) = cj;
    FhRhd = reshape(FhRhd, [Nkx Nky Nkz Nc]);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %% Calculate B^H * S^H * (I_{Nc} kron F^H * R_{Omega}^H) * d
    outp = sum(bsxfun(@times, conj(sens(:,:,:,:,nacq)), FhRhd), 4); % Nkx x Nky x Nkz

    %% Bochao: Add phase constraints
    outp = bsxfun(@times, outp, exp(-1i*phase(:,:,nacq)));

    %% Vectorize the output
    outp = reshape(outp, [N 1]); % Nk x Nky x Nkz => N x 1

    %% Tikhonov regularization
    outp = outp + lambda * inp;
    
    outp_dual = outp_dual + outp;
end

end