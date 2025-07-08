function outp_dual = type1_pec_sense_adjoint_dual(inp, sens, mask, nj, p1, p2, iflag, eps, support_mask, topup_sign, phase)

%% Start a stopwatch timer
start_time = tic;

%% Get imaging parameters
[Nkx,Nky,Nkz,Nc,Nacq] = size(sens);

%--------------------------------------------------------------------------
% Calculate the total number of voxels
%--------------------------------------------------------------------------
N = Nkx * Nky * Nkz;

%% Calculate scale factors
type2_scale_factor = 1 / sqrt(Nkx * Nky);

outp_dual = zeros(N,1);

for nacq = 1:Nacq
    %% Calculate (I_{Nc} kron R_{Omega}^H) * d
    Rhd = bsxfun(@times, mask(:,:,nacq), inp(:,:,:,:,nacq));

    %% Calculate (I_{Nc} kron F^H) * (I_{Nc} kron R_{Omega}^H) * d
    %--------------------------------------------------------------------------
    % Perform type-2 NUFFT (nonuniform, image space <= uniform, k-space)
    % Siemens: k-space <=> image space
    %--------------------------------------------------------------------------
    tstart = tic; fprintf('%s: Performing type-2 NUFFT... ', datetime);
    FhRhd = complex(zeros(N, Nc, 'single'));
    cj = type2_scale_factor * finufft2d2(p1{nacq}, p2{nacq}, iflag, eps, reshape(double(Rhd(:,:,1,:)), [Nkx Nky Nc])); % N x Nc
    FhRhd((support_mask(:,:,nacq) > 0),:) = cj;
    FhRhd = reshape(FhRhd, [Nkx Nky Nkz Nc]);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %% Calculate B^H * S^H * (I_{Nc} kron F^H * R_{Omega}^H) * d
    outp = sum(bsxfun(@times, conj(sens(:,:,:,:,nacq)), FhRhd), 4); % [Nkx x Nky x Nkz]

    %% Bochao: Add phase constraints
    outp = bsxfun(@times, outp, exp(-1i*phase(:,:,nacq)));

    %% Vectorize the output
    outp = reshape(outp, [N 1]); % Nk x Nky x Nkz => N x 1

    outp_dual = outp_dual + outp;
end
end