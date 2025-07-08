function output = gradstep(I1,I2,Tx,Btotal,Bgrad)

[~,gradTx,~] = gradient(Tx);

I1mod = movepixels_3d(I1,-Tx);
I2mod = movepixels_3d(I2,Tx);

%gradient of image.
I1grad = imagegrad(I1,-Tx);
I2grad = imagegrad(I2,Tx);

imagediff = (1 - gradTx) .* I1mod - (1 + gradTx) .* I2mod;

im1 = imagediff .* I1grad .* (1 - gradTx);
im2 = imagediff .* I1mod;
im3 = imagediff .* I2grad .* (1 + gradTx);
im4 = imagediff .* I2mod;

df1 = Btotal * (im1(:) + im3(:));
df2 = Bgrad * (im2(:) + im4(:));

output = -(df1 + df2);


