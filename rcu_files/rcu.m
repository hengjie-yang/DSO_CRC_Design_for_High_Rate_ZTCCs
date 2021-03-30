function [Pe,rhop] = rcu(n, R, pX, X, sigma2)


% briefly check, where we are operating
[~,Rcrit]=E0funs(pX, X, sigma2, 1);
Rcrit = Rcrit/log(2);

[~,Rmi]=E0funs(pX, X, sigma2, 0);
Rmi = Rmi/log(2);

% operating between Rcrit and MI
if (R >= Rcrit) && (R <= Rmi)
%    fprintf('SNR: %.2f: Operating within Rcrit = %.2f and MI = %.2f.\n', 10*log10(1/sigma2), Rcrit, Rmi);
   mode = 0;
   rhop = fzero(@(x) find_rho(x), [0 1]);
elseif R <= Rcrit
%    fprintf('SNR: %.2f: Operating below Rcrit = %.2f (MI = %.2f).\n', 10*log10(1/sigma2), Rcrit, Rmi);
   mode = 1;
   rhop = fzero(@(x) find_rho(x), [1 100]);
else
%    fprintf('SNR: %.2f: Operating above MI = %.2f (Rcrit = %.2f).\n', 10*log10(1/sigma2), Rmi, Rcrit);
   mode = 2;
   rhop = fzero(@(x) find_rho(x), [-0.98 0.1]); % original [-1.99 0.1]
end


[E0,E0p,E0pp] = E0funs(pX,X,sigma2,rhop);
V = -E0pp;

% omegapp = integral(@(y) Qpmf(y,pX,X,sigma2,rhop) .* omegapp_fun(y,rhop), -20, 20)
omegapp = integral(@(y) fun(y,rhop), -20, 20);


if mode == 0 % medium SNR region
    
    theta_n = 1/sqrt(1+rhop) * ((1+rhop)/(sqrt(2*pi*n*omegapp)))^rhop;
    Pe = theta_n * (Psi(rhop*sqrt(n*V)) + Psi((1-rhop) * sqrt(n*V))) * 2^(-n*(E0_integral(pX,X,sigma2,rhop)/log(2) - rhop * R));
elseif mode == 1 % high SNR region
    
    theta_n = 1/sqrt(1+rhop) * ((1+rhop)/(sqrt(2*pi*n*omegapp)))^rhop;
    omegapp1 = integral(@(y) fun(y,1), -20, 20);
    theta_n1 = 1/sqrt(1+1) * ((1+1)/(sqrt(2*pi*n*omegapp1)))^1;
    [E01] = E0funs(pX,X,sigma2,1);
    Pe =  exp(-n*(E01-R*log(2))) * theta_n1 + theta_n * (Psi(rhop*sqrt(n*V)) + Psi((1-rhop) * sqrt(n*V))) * 2^(-n*(E0_integral(pX,X,sigma2,rhop)/log(2) - rhop * R));
else
    
    theta_n = 1/sqrt(1+rhop) * ((1+rhop)/(sqrt(2*pi*n*omegapp)))^rhop;
    Pe = 1 + theta_n * (Psi(rhop*sqrt(n*V)) + Psi((1-rhop) * sqrt(n*V))) * 2^(-n*(E0_integral(pX,X,sigma2,rhop)/log(2) - rhop * R));
end

    function val = fun(y,rho)
        val = Qpmf(y,pX,X,sigma2,rho) .* omegapp_fun(y,rho);
    end

    function val = find_rho(x)
        val = E0prime_integral(pX,X,sigma2,x)/log(2) - R;
    end

    function val = Psi(x)
        val = .5*erfc(abs(x)/sqrt(2))*exp(x^2/2)*sign(x);
    end

    function val = omegapp_fun(y,rho)
        tau = 1/(1+rho);
        y_s = size(y);
        y = y(:);
        Xtmp = reshape(X,1,length(X));
        pXtmp = reshape(pX,1,length(X));
        val = (sum(normpdf(y,Xtmp,sqrt(sigma2)).^tau.*pXtmp, 2) .* sum(normpdf(y,Xtmp,sqrt(sigma2)).^tau.*pXtmp.*log(normpdf(y,Xtmp,sqrt(sigma2))).^2, 2) - sum(normpdf(y,Xtmp,sqrt(sigma2)).^tau.*pXtmp.*log(normpdf(y,Xtmp,sqrt(sigma2))), 2).^2) ./ ...
           sum(normpdf(y,Xtmp,sqrt(sigma2)).^tau.*pXtmp, 2).^2;
        val = reshape(val,y_s);
        val(isnan(val)) = 0;
    end


end


