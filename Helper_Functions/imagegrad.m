function output = imagegrad(Iin,Tx)

[~,gradX,~] = gradient(Iin);

gradX = movepixels_3d(gradX,Tx);

output = gradX;