function output = bsplines(n,h)

%n is the grid size
%h is the spline width (essentially). a spline covers 4h + 1 knots.
%m is the number of splines

x = 1:n;
knotcenters = 1:h:n;

output = zeros(n,length(knotcenters));
for index = -ceil(2*h):1:ceil(2*h)
    dist = abs(knotcenters - x.')/h;
    %imagesc(dist);
        
    kk = abs(dist) <= 1;
    output(kk) = 2/3 - (1 - dist(kk)/2).*dist(kk).^2;
    
    kk = (abs(dist) > 1 & abs(dist) <= 2);
    output(kk) = ((2 - dist(kk)).^3)/6;
end