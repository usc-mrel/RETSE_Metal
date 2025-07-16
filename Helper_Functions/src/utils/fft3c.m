function res = fft3c(x)

S = size(x);
fctr = S(1)*S(2)*S(3);

x = reshape(x,S(1),S(2),S(3),prod(S(4:end)));

res = zeros(size(x));
for n=1:size(x,4)
	res(:,:,:,n) = 1/sqrt(fctr)*fftshift(fftn(ifftshift(x(:,:,:,n))));
end

res = reshape(res,S);



