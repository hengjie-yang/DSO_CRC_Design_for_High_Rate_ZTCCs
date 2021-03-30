function val = E0_integral(pX,X,sigma2,rho)

val = -log(integral(@(y) inner_func(y), -20, 20));


    function val1 = inner_func(yvals)
        val1 = 0;
        yvals_s = size(yvals);
        yvals = yvals(:);
        for ii=1:length(X)
            val1 = val1 + 1/sqrt(2*pi*sigma2) * exp(-(yvals-X(ii)).^2/(2*sigma2)) * pX(ii) .* exp(-rho * inf_density(X(ii), yvals, pX, X, sigma2, 1/(1+rho)));
        end
        val1 = reshape(val1,yvals_s);
        val1(isnan(val1)) = 0;
    end

end
