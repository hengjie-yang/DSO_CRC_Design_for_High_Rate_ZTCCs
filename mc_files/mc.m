function [Pe, rhop] = mc(n, R, pX, X, sigma2)

%   The function is to evaluate the saddlepoint approximation of the
%   meta-converse bound in [1], Eq. (42) for the BI-AWGN channel.
%
%   Reference:
%   [1] J. Font-Segura et al., "Saddlepoint approximations of lower and
%       upper bounds to the error probability in channel coding."
%
%   Input parameters:
%       1) n: a scalar denoting the blocklength
%       2) R: a scalar denoting the channel coding rate
%       3) pX: a column vector denoting the unif. distribution on constellation
%       symbols
%       4) X: a column vector denooting the constellation symbols
%       5) sigma2: a scalar denoting the AWGN noise variance.
%
%   Output parameters:
%       1) Pe: the meta-converse value given by saddlepoint approximation.
%       2) rhop: the optimal parameter \rho that yields Pe.
%
%   Written by Hengjie Yang (hengjie.yang@ucla.edu)     01/05/21.
%


% briefly check, where we are operating
[~,Rcrit]=E0funs(pX, X, sigma2, 1);
Rcrit = Rcrit/log(2);

[~,Rmi]=E0funs(pX, X, sigma2, 0);
Rmi = Rmi/log(2);

if (R >= Rcrit) && (R <= Rmi)
   fprintf('SNR: %.2f: Operating within Rcrit = %.2f and MI = %.2f.\n', 10*log10(1/sigma2), Rcrit, Rmi);
elseif R <= Rcrit
   fprintf('SNR: %.2f: Operating below Rcrit = %.2f (MI = %.2f).\n', 10*log10(1/sigma2), Rcrit, Rmi);
else
   fprintf('SNR: %.2f: Operating above MI = %.2f (Rcrit = %.2f).\n', 10*log10(1/sigma2), Rmi, Rcrit);
end



[rhop, Pe] = fminbnd(@fun, 0, 2);
Pe = -Pe;

    function val = fun(rho)
        [E0, E0p, E0pp] = E0funs(pX, X, sigma2, rho);       
        U = -(1+rho)*E0pp; % after simplifying Eq. (43) in Reference [1].
        val = -exp(-n*(E0 - rho*E0p))*(Psi(sqrt(n*U))+Psi(rho*sqrt(n*U))-exp(-n*(R*log(2)-E0p)));
    end
    
    function val = Psi(x)
        val = .5*erfc(abs(x)/sqrt(2))*exp(x^2/2)*sign(x);
    end


% rhos = 0:0.1:10;
% vals = zeros(1, size(rhos, 2));
% for ii = 1:length(rhos)
%     vals(ii) = fun(rhos(ii));
% end
% plot(rhos, vals, '--');
% grid on


end







