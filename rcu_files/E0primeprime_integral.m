function val = E0prime_integral(pX,X,sigma2,rho)

X = reshape(X,1,length(X));
pX = reshape(pX,1,length(pX));

num = integral(@(y) inner_func(y), -30, 30);
  

    function val1 = inner_func(yvals)
        val1 = 0;
        yvals_s = size(yvals);
        yvals = yvals(:);
        for ii=1:length(X)
            val1 = val1 + 1/sqrt(2*pi*sigma2) * exp(-(yvals-X(ii)).^2/(2*sigma2)) * pX(ii) .* ...
                exp(-rho * inf_density(X(ii), yvals, pX, X, sigma2, 1/(1+rho))) .* ...
                inf_density(X(ii), yvals, pX, X, sigma2, 1/(1+rho)); % .* ...
               % (sum(normpdf(yvals,X,sqrt(sigma2)).^(1/(1+rho)).*pX,2).*log(normpdf(yvals,X(ii),sqrt(sigma2))) - sum(normpdf(yvals,X,sqrt(sigma2)).^(1/(1+rho)).*pX.*log(normpdf(yvals,X,sqrt(sigma2))),2) ./ sum(normpdf(yvals,X,sqrt(sigma2)).^(1/(1+rho)) .* pX,2)) * (-1*(1+rho))^(-2);
        end
        val1 = reshape(val1,yvals_s);
        val1(isnan(val1)) = 0;
    end


denom = exp(-E0_integral(pX,X,sigma2,rho));

val = num/denom;

end