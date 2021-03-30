function res = inf_density(x,y, pX, X, sigma2, s)

y = y(:);
% x = x(:);

X = reshape(X, 1,length(X));
pX = reshape(pX, 1,length(X));

num = exp(-s*(y-x).^2/(2*sigma2));
denom = sum(exp(-s*(y-X).^2/(2*sigma2)) .* pX, 2);

res = log(num./denom);
 % catch cases where we have sth like log(0/0); this may happen if we
 % evaluate the density in a y which is far from probable. however, this
 % might happen depending on the context, i.e., when integral makes a
 % particular bad choice for the intermediate points
res(isnan(res)) = 0;
% catch cases where we have a very small numerator, i.e., s.th. like
% log(0/epsilon)
res(res==-inf) = -999;
% res(res==+inf) = +999;

end
