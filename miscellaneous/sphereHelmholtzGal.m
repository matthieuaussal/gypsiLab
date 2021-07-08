function u = sphereHelmholtzGal(name,bndCond,rho,k,X,type,X0) 
% Copyright (c) 20015-2017, Matthieu Aussal, Ecole Polytechnique       
% GNU General Public License v3.0. 
% Computation of analytic field for helmoltz spherical scattering

% Sign convention
if k <= 0
    sgnk = @(v) v;
else
    sgnk = @(v) conj(v);
end 
k = abs(k);

% Spherical data
[phi,theta,r] = cart2sph(X(:,1),X(:,2),X(:,3));
phi = phi'; theta = pi/2 - theta';  r = r';
[phi_inc,theta_inc,~] = cart2sph(X0(:,1),X0(:,2),X0(:,3));
theta_inc = pi/2 - theta_inc;
kR = k.*r;
u = zeros(1,length(theta));

% Boundary condition
if strcmp(bndCond,'dir') % Dirichlet
    jsph = @(n,z) sqrt(pi./(2.*z)) .* besselj(n+0.5,z); % Bessel
    hsph = @(n,z) sqrt(pi./(2.*z)) .* besselh(n+0.5,z); % Hankel
    alpha = @(n,z) - jsph(n,z)./hsph(n,z); % Scattering coefficient
elseif strcmp(bndCond,'neu') % Neumann
    jsph = @(n,z) sqrt(pi./(2.*z)) .* besselj(n+0.5,z);
    hsph = @(n,z) sqrt(pi./(2.*z)) .* besselh(n+0.5,z);
%     djsph = @(n,z) (n.*jsph(n-1,z) - (n+1).*jsph(n+1,z))/(2*n+1) ;    
%     dhsph = @(n,z) -(1i./z.^2 - djsph(n,z).*hsph(n,z))./jsph(n,z);
    djsph = @(n,z) n./z.*jsph(n,z) - jsph(n+1,z);    
    dhsph = @(n,z) n./z.*hsph(n,z) - hsph(n+1,z);
    alpha = @(n,z) - djsph(n,z)./dhsph(n,z);
end

% Harmonic coefficient
if strcmp(type,"plane")
    coeff = @(n) 4*pi*(1i)^n;
elseif strcmp(type,"spher")
    coeff = @(n) 4*pi*1i*k*hsph(n,k*norm(X0));
end

% Infinite spherical radiation
if strcmp(name,'inf')   
    n = 0;
    while (max(abs(alpha(n,kR))) > 1e-12)
        Pn = legendre(n,cos(theta),"norm");
        Pn_inc = legendre(n,cos(theta_inc),"norm");
        for m=0:n
            % Associated Legendre polynomial of degree n and order m
            Pnm = Pn(m+1,:);       % For the whole domain
            Pnm_inc = Pn_inc(m+1); % For the incident wave
            
            % Spherical Harmonic of degree n and order m
            a = 2*n+1;
            b = 4*pi*(n+1/2);
            C = ((-1)^m)*sqrt(a/b);
            Ynm = C .* Pnm .* exp(1i*m*phi);
            Ynm_inc = C .* Pnm_inc .* exp(1i*m*phi_inc);

            % Add coefficient to sum
            u = u + (-1i)^(n)*coeff(n).*alpha(n,k*rho)...
              .* ((m==0)*Ynm.*conj(Ynm_inc)+2.*(m>0).*real(conj(Ynm).*Ynm_inc));
        end
        n = n+1;
    end
    u = -conj((1i/k) .* u.');
    
% Finite spherical radiation
elseif strcmp(name,'bnd') || strcmp(name,'dom')
    n = 0;
    hn = hsph(n+1,kR);
    while (max(abs(alpha(n,kR).*hn)) > 1e-12) + (n<10)
        Pn = legendre(n,cos(theta),"norm");
        Pn_inc = legendre(n,cos(theta_inc),"norm");
        hn = hsph(n,kR);
        for m=0:n
            % Associated Legendre polynomial of degree n and order m
            Pnm = Pn(m+1,:);       % For the whole domain
            Pnm_inc = Pn_inc(m+1); % For the incident wave
            
            % Spherical Harmonic of degree n and order m
            a = 2*n+1;
            b = 4*pi*(n+1/2);
            C = ((-1)^m)*sqrt(a/b);
            Ynm = C .* Pnm .* exp(1i*m*phi);
            Ynm_inc = C .* Pnm_inc .* exp(1i*m*phi_inc);

            % Add coefficient to sum
            u = u + coeff(n).*(alpha(n,k*rho).*hn)...
              .* ((m==0)*Ynm.*conj(Ynm_inc)+2.*(m>0).*real(conj(Ynm).*Ynm_inc));
        end
        n = n+1;
    end
    u = conj(u.');
    if strcmp(name,'dom')
        u(r<rho) = 0;
    end
    
else
    error('sphereHelmholtz.m : unavailable case')
end

% output
u = sgnk(u);
end