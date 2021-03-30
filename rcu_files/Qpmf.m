function q = Qpmf(y,pX,X,sigma2,rho)

X = reshape(X,1,length(X));
pX = reshape(pX,1,length(pX));

C = integral(@(yvals) fun(yvals), -20, 20);

q = fun(y) / C;

    function val1 = fun(yvals)
        yvals_s = size(yvals);
        yvals = yvals(:);
%         for ii=1:length(X)
%             val1 = val1 + normpdf(yvals,X(ii),sqrt(sigma2)).*pX(ii).*exp(-rho*inf_density(X(ii),yvals,pX,X,sigma2,1/(1+rho)));
%         end
        val1 = sum(normpdf(yvals,X,sqrt(sigma2)).*pX.*exp(-rho*inf_density(X,yvals,pX,X,sigma2,1/(1+rho))),2);
        val1(isnan(val1)) = 0;
        val1 = reshape(val1,yvals_s);
    end

end