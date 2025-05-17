% This code accompanies the paper:
%      Wu, B., 2024. An Extension of the Euler-Maclaurin Summation Formula 
%      to Nearly Singular Functions. arXiv preprint arXiv:2409.19192.
%
% Four examples are presented in this code:
% Examples 1 & 2 use analytic expressions for near-singular correction. (cf. the arXiv paper)
% Examples 1b & 2b use finite difference approximation for truncated near-singular correction.

% NOTE: Symbolic Toolbox is required to run this code.

%% Example 1: Integrate f(x) = d*exp(x)/(c^2*x^2 + d^2) from -1 to 1 using trapezoidal rule
% NOTE: Correction computed using analytic expression
% Exact solution: anti-derivative is 
%                    imag(exp(-1i*d/c)/c*[expint(-x-1i*d/c)-expint(x-1i*d/c)])
%                 or imag(exp(-1i*d/c)/c*[Ei(-x+1i*d/c)-Ei(x+1i*d/c)]
% WARNING: Matlab's expint(x) is E_1(x)==-Ei(-x), not Ei(x).

dvals = [0.1, 0.01, 0.0001];
c = 1;

for kplot = 1:3
    d = dvals(kplot);
    Iexact = imag( exp(-1i*d/c)/c * (expint(-1-1i*d/c) - expint(1-1i*d/c)) );

    err_nat = [];
    err_corr = [];
    N = 10*2.^(1:0.25:10); N = round(N/2)*2;
    for n = N

        % punctured trapezoidal rule
        x = linspace(-1,1,n+1);
        h = 2/n;
        w = h*ones(size(x));
        w0 = w;
        w0(n/2+1) = 0;

        % edge corrections via Gregory quadrature
        p = 12; g = gregory_weights(p);
        w(1:p-1) = w(1:p-1) + g*h;
        w(end:-1:end-p+2) = w(end:-1:end-p+2) + g*h;
        w0(1:p-1) = w0(1:p-1) + g*h;
        w0(end:-1:end-p+2) = w0(end:-1:end-p+2) + g*h;

        % evaluate near-singular integral on [-1,1]
        f = @(x) d*exp(x)./(c^2*x.^2+d^2);
        f_top = @(x) d*exp(x);  % smooth numerator of integrand
        I_nat = dot(w,f(x));    % trapezoidal rule
        I0 = dot(w0,f(x));      % punctured trapezoidal rule

        % near-singular correction (using analytic expression)
        lam = d/c/h;
        z0 = imag(double(psi(sym(1)+1i*lam)))/lam;
        E0 = 1/c^2*real(f_top(1i*lam*h)-f_top(0))/(-lam^2*h^2)*h + (pi/c/d - 2*z0/c^2/h) * real(f_top(1i*lam*h));
        Icorr = I0 + E0;

        err_nat = [err_nat, abs(I_nat - Iexact)];
        err_corr = [err_corr, abs(Icorr - Iexact)];
    end

    figure(1)
    subplot(2,3,kplot)
    plot(x,f(x),'LineWidth',2)
    title(['$d = ',num2str(d,'%g'),'$'],'Interpreter','latex')
    if kplot == 1
    ylabel('$f(x)$','Interpreter','latex','Rotation',0)
    end

    subplot(2,3,kplot+3)
    loglog(N,err_nat,'o-'), hold on
    loglog(N,err_corr,'*-'), hold off
    ylim([1e-16,1e3])
    yticks([1e-15,1e-8,1e0])
    xlabel('$N$','Interpreter','latex')
    if kplot == 1
        ylabel('error','Interpreter','latex')
        legend('no correction','corrected')
    end
end

%% Example 2: Integrate f(x) = d*exp(x)/(c^2(x-xs)^2 + d^2) from -1 to 1 using trapezoidal rule
% NOTE: Correction computed using analytic expression
% Exact solution: anti-derivative is 
%                    imag(exp(s-1i*d/c)/c*[expint(-x+s-1i*d/c)-expint(x+s-1i*d/c)])
%                 or imag(exp(s-1i*d/c)/c*[Ei(-x-s+1i*d/c)-Ei(x-s+1i*d/c)]
% WARNING: Matlab's expint(x) is E_1(x)==-Ei(-x), not Ei(x).

dvals = [0.1, 0.01, 0.0001];
xs = 0.1;
c = 1.21;

for kplot = 1:3
    d = dvals(kplot);
    Iexact = imag( exp(xs-1i*d/c)/c * (expint(-1+xs-1i*d/c) - expint(1+xs-1i*d/c)) );
    
    err_nat = [];
    err_corr = [];
    N = 10*2.^(1:0.25:10); N = round(N/2)*2;
    for n = N

        % punctured trapezoidal rule
        x = linspace(-1,1,n+1);
        h = 2/n;
        w = h*ones(size(x));
        ind = find(abs(x-xs)==min(abs(x-xs)),1,"first"); % closest point
        w0 = w;
        w0(ind) = 0;
        x0 = x(ind);
        s0 = (xs-x0)/h;

        % edge corrections via Gregory quadrature
        p = 12; g = gregory_weights(p);
        w(1:p-1) = w(1:p-1) + g*h;
        w(end:-1:end-p+2) = w(end:-1:end-p+2) + g*h;
        w0(1:p-1) = w0(1:p-1) + g*h;
        w0(end:-1:end-p+2) = w0(end:-1:end-p+2) + g*h;

        % evaluate near-singular integral on [-1,1]
        f = @(x) d*exp(x)./(c^2*(x-xs).^2+d^2);
        f_top = @(x) d*exp(x);  % smooth numerator of integrand
        I_nat = dot(w,f(x));    % trapezoidal rule
        I0 = dot(w0,f(x));      % punctured trapezoidal rule

        % near-singular correction (using analytic expression)
        lam = d/c/h;
        p0s = -imag(double(psi(sym(1)-s0-1i*lam)+psi(sym(1)+s0-1i*lam)))/lam;
        p1s = -real(double(psi(sym(1)-s0-1i*lam)-psi(sym(1)+s0-1i*lam)));
        E0 =  -1/(h*c^2) * ( - f_top(x0)/(lam^2+s0^2) ...
                             + (p0s+ 1/(s0^2+lam^2) - pi/lam) * real(f_top(xs+1i*lam*h)) ...
                             + (p1s-s0/(s0^2+lam^2)) * imag(f_top(xs+1i*lam*h)/lam) );
        Icorr = I0 + E0;
        
        err_nat = [err_nat, abs(I_nat - Iexact)];
        err_corr = [err_corr, abs(Icorr - Iexact)];
    end

    figure(1)
    subplot(2,3,kplot)
    plot(x,f(x),'LineWidth',2)
    title(['$d = ',num2str(d,'%g'),'$'],'Interpreter','latex')
    if kplot == 1
    ylabel('$f(x)$','Interpreter','latex','Rotation',0)
    end

    subplot(2,3,kplot+3)
    loglog(N,err_nat,'o-','MarkerSize',5), hold on
    loglog(N,err_corr,'*-','MarkerSize',5), hold off
    ylim([1e-16,1e3])
    yticks([1e-15,1e-8,1e0])
    xlabel('$N$','Interpreter','latex')
    if kplot == 1
        ylabel('error','Interpreter','latex')
        legend('no correction','corrected')
    end
end

%% Example 1b: Integrate f(x) = d*exp(x)/(c^2 x^2 + d^2) from -1 to 1 using trapezoidal rule
% NOTE: Correction weights computed using finite difference
% Exact solution: anti-derivative is 
%                    imag(exp(-1i*d/c)/c*[expint(-x-1i*d/c)-expint(x-1i*d/c)])
%                 or imag(exp(-1i*d/c)/c*[Ei(-x+1i*d/c)-Ei(x+1i*d/c)]
% WARNING: Matlab's expint(x) is E_1(x)==-Ei(-x), not Ei(x).

dvals = [0.1, 0.01, 0.0001];
c = 1.21;

for kplot = 1:3
    d = dvals(kplot);
    Iexact = imag( exp(-1i*d/c)/c * (expint(-1-1i*d/c) - expint(1-1i*d/c)) );

    err_nat = [];
    err_corr = [];
    N = 10*2.^(1:0.25:10); N = round(N/2)*2;
    for n = N

        % punctured trapezoidal rule
        x = linspace(-1,1,n+1);
        h = 2/n;
        w = h*ones(size(x));
        w0 = w;
        w0(n/2+1) = 0;

        % edge corrections via Gregory quadrature
        p = 12; g = gregory_weights(p);
        w(1:p-1) = w(1:p-1) + g*h;
        w(end:-1:end-p+2) = w(end:-1:end-p+2) + g*h;
        w0(1:p-1) = w0(1:p-1) + g*h;
        w0(end:-1:end-p+2) = w0(end:-1:end-p+2) + g*h;

        % evaluate near-singular integral on [-1,1]
        f = @(x) d*exp(x)./(c^2*x.^2+d^2);
        f_top = @(x) d*exp(x);  % smooth numerator of integrand
        I_nat = dot(w,f(x));    % trapezoidal rule
        I0 = dot(w0,f(x));      % punctured trapezoidal rule

        % near-singular correction (using finite difference)
        m = 3; % O(h^2m) correction
        lam = d/c/h;
        w = digamma_weights_OnMesh(m,lam); % weights on (half) stencil 0:m
        w = [flipud(w(2:end));w]/(c^2*h);  % weights on full stencil -m:m by symmetry
        E0 = f_top((-m:m)*h)*w;
        Icorr = I0 + E0;

        err_nat = [err_nat, abs(I_nat - Iexact)];
        err_corr = [err_corr, abs(Icorr - Iexact)];
    end

    figure(1)
    subplot(2,3,kplot)
    plot(x,f(x),'LineWidth',2)
    title(['$d = ',num2str(d,'%g'),'$'],'Interpreter','latex')
    if kplot == 1
    ylabel('$f(x)$','Interpreter','latex','Rotation',0)
    end

    subplot(2,3,kplot+3)
    loglog(N,err_nat,'o-'), hold on
    loglog(N,err_corr,'*-'), hold off
    ylim([1e-16,1e3])
    yticks([1e-15,1e-8,1e0])
    xlabel('$N$','Interpreter','latex')
    if kplot == 1
        ylabel('error','Interpreter','latex')
        legend('no correction','corrected')
    end
end

%% Example 2b: Integrate f(x) = d*exp(x)/(c^2(x-xs)^2 + d^2) from -1 to 1 using trapezoidal rule
% NOTE: Correction weights computed using finite difference
% Exact solution: anti-derivative is 
%                    imag(exp(s-1i*d/c)/c*[expint(-x+s-1i*d/c)-expint(x+s-1i*d/c)])
%                 or imag(exp(s-1i*d/c)/c*[Ei(-x-s+1i*d/c)-Ei(x-s+1i*d/c)]
% WARNING: Matlab's expint(x) is E_1(x)==-Ei(-x), not Ei(x).

dvals = [0.1, 0.01, 0.0001];
xs = 0.1;
c = 1.21;

for kplot = 1:3
    d = dvals(kplot);
    Iexact = imag( exp(xs-1i*d/c)/c * (expint(-1+xs-1i*d/c) - expint(1+xs-1i*d/c)) );
    
    err_nat = [];
    err_corr = [];
    N = 10*2.^(1:0.25:10); N = round(N/2)*2;
    for n = N

        % punctured trapezoidal rule
        x = linspace(-1,1,n+1);
        h = 2/n;
        w = h*ones(size(x));
        ind = find(abs(x-xs)==min(abs(x-xs)),1,"first"); % closest point
        w0 = w;
        w0(ind) = 0;
        x0 = x(ind);
        s0 = (xs-x0)/h;

        % edge corrections via Gregory quadrature
        p = 12; g = gregory_weights(p);
        w(1:p-1) = w(1:p-1) + g*h;
        w(end:-1:end-p+2) = w(end:-1:end-p+2) + g*h;
        w0(1:p-1) = w0(1:p-1) + g*h;
        w0(end:-1:end-p+2) = w0(end:-1:end-p+2) + g*h;

        % evaluate near-singular integral on [-1,1]
        f = @(x) d*exp(x)./(c^2*(x-xs).^2+d^2);
        f_top = @(x) d*exp(x);  % smooth numerator of integrand
        I_nat = dot(w,f(x));    % native quadrature
        I0 = dot(w0,f(x));      % punctured trapezoid

        % near-singular correction (using finite difference)
        m = 3; % O(h^2m) correction
        lam = d/c/h;
        w = digamma_weights_Offmesh(m,lam,s0)/(c^2*h); % weights on stencil -m:m
        E0 = f_top(x0+(-m:m)*h)*w;
        Icorr = I0 + E0;
        
        err_nat = [err_nat, abs(I_nat - Iexact)];
        err_corr = [err_corr, abs(Icorr - Iexact)];
    end

    figure(1)
    subplot(2,3,kplot)
    plot(x,f(x),'LineWidth',2)
    title(['$d = ',num2str(d,'%g'),'$'],'Interpreter','latex')
    if kplot == 1
    ylabel('$f(x)$','Interpreter','latex','Rotation',0)
    end

    subplot(2,3,kplot+3)
    loglog(N,err_nat,'o-','MarkerSize',5), hold on
    loglog(N,err_corr,'*-','MarkerSize',5), hold off
    ylim([1e-16,1e3])
    yticks([1e-15,1e-8,1e0])
    xlabel('$N$','Interpreter','latex')
    if kplot == 1
        ylabel('error','Interpreter','latex')
        legend('no correction','corrected')
    end
end


function w = digamma_weights_OnMesh(m,lam)
% Near-singular correction weights using centered difference formula.
% Near singularity is on mesh.
% Input
%   m : to compute O(h^2m) correction weights
%   lam : parameter for correction 
% Output
%   w : correction weights on the (half) stencil 0:m
%   NOTE: weights on the remaining stencil -m:-1 can be obtained using even symmetry

z = -imag(double(psi(sym(1)-1i*lam)))/lam;
l = (-lam^2).^(0:m)';
b = l*(-2*z + pi/lam) + [0;l(1:m)];

A = 2*(0:m).^((0:2:2*m)'); A(1) = 1;
w = A\b;
end

function w = digamma_weights_Offmesh(m,lam,s)
% Near-singular correction weights using centered difference formula.
% Near singularity is off mesh.
% Input
%   m : to compute O(h^2m) correction weights
%   lam : parameter for correction 
%   s : parameter for off-mesh shift
% Output
%   w : correction weights on the stencil -m:m

p(1) = -imag(double(psi(sym(1)-s-1i*lam)+psi(sym(1)+s-1i*lam)))/lam;
if m > 0
    p(2) = -real(double(psi(sym(1)-s-1i*lam)-psi(sym(1)+s-1i*lam)));
    for k = 0:2*m-2
        p(k+3) = -(-s)^k-lam^2*p(k+1);
    end
end
b = -p;
b(1:2:end) = b(1:2:end) + (-1).^(0:m).*lam.^(-1:2:2*m-1)*pi;

b = b(:);
A = ((-m:m)-s).^((0:2*m)');
w = A\b;
end

function d=gregory_weights(p)
% Edge correction weights for the trapezoidal rule with weight h=1
% Uncorrected weights = all 1's
% Corrected weights = 1+d
% p = order

switch p
    case 2
        d = -1/2;
    case 3
        d = [-7/12, 1/12];
    case 4
        d = [-5/8, 1/6, -1/24];
    case 5
        d = [-469/720, 59/240, -29/240, 19/720];
    case 6
        d = [-193/288, 77/240, -7/30, 73/720, -3/160];
    case 7
        d = [-41393/60480, 23719/60480, -11371/30240, 7381/30240, -5449/60480, 863/60480];
    case 8
        d = [-12023/17280, 6961/15120, -66109/120960, 33/70, -31523/120960, 1247/15120, -275/24192];
end

% note: when p >= 10, some of the modified weights (1+d) are negative,
% quadrature may become unstable
if p > 8
    M = p-2;
    P = (0:M)';
    rhs = -zeta(-P);
    rhs(1) = -rhs(1);
    A = vpa(0:M).^P;
    d = double(A\rhs)';
end
end
