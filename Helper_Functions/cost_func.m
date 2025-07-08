function output = cost_func(I1,I2,Tx)

[~,gradTx,~] = gradient(Tx);

I1mod = (1 - gradTx) .* movepixels_3d(I1,-Tx);
I2mod = (1 + gradTx) .* movepixels_3d(I2,Tx);

output = sum(abs(I1mod(:) - I2mod(:)).^2); 