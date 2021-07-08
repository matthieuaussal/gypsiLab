%+========================================================================+
%|                                                                        |
%|            This script uses the GYPSILAB toolbox for Matlab            |
%|                                                                        |
%| COPYRIGHT : Matthieu Aussal (c) 2015-2017.                             |
%| PROPERTY  : Centre de Mathematiques Appliquees, Ecole polytechnique,   |
%| route de Saclay, 91128 Palaiseau, France. All rights reserved.         |
%| LICENCE   : This program is free software, distributed in the hope that|
%| it will be useful, but WITHOUT ANY WARRANTY. Natively, you can use,    |
%| redistribute and/or modify it under the terms of the GNU General Public|
%| License, as published by the Free Software Foundation (version 3 or    |
%| later,  http://www.gnu.org/licenses). For private use, dual licencing  |
%| is available, please contact us to activate a "pay for remove" option. |
%| CONTACT   : matthieu.aussal@polytechnique.edu                          |
%| WEBSITE   : www.cmap.polytechnique.fr/~aussal/gypsilab                 |
%|                                                                        |
%| Please acknowledge the gypsilab toolbox in programs or publications in |
%| which you use it.                                                      |
%|________________________________________________________________________|
%|   '&`   |                                                              |
%|    #    |   FILE       : nrtHmxHelmholtzBWdir.m                        |
%|    #    |   VERSION    : 0.40                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 14.03.2018                                    |
%| ( === ) |   SYNOPSIS   : Solve dirichlet scatering problem with        |
%|  `---'  |                Brackage-Werner formulation                   |
%+========================================================================+

% Cleaning
clear all
close all
clc

% Gypsilab path
run('../../addpathGypsilab.m')

% Parameters
N   = 1e3
tol = 1e-3
typ = 'P1'
gss = 3
X0  = [0 0 -1];
type = "plane"; % "plane" / "spher"

% Spherical mesh
sphere = mshSphere(N,1);
sigma  = dom(sphere,gss);    
figure
plot(sphere)
axis equal

% Radiative mesh
square     = mshSquare(5*N,[5 5]);
square.vtx = [square.vtx(:,1) zeros(size(square.vtx,1),1) square.vtx(:,2)];
hold on
plot(square)

% Frequency adjusted to maximum esge size
stp = sphere.stp;
k   = 1/stp(2)
f   = (k*340)/(2*pi)

% Incident wave
PW = @(X) exp(1i*k*X*(X0'./norm(X0)));
Norme = @(X1,X2) sqrt( (X1(:,1)-X2(1)).^2 + (X1(:,2)-X2(2)).^2 + (X1(:,3)-X2(3)).^2 );
SW = @(X) exp(1i*k*Norme(X,X0))./Norme(X,X0);
if strcmp(type,"plane")
    IW = @(X) PW(X);
elseif strcmp(type,"spher")
    IW = @(X) SW(X);
else
    error('nrtHelmholtzGalAnalyticPWdir.m : unavailable case')
end

% Incident wave representation
plot(sphere,real(IW(sphere.vtx)))
plot(square,real(IW(square.vtx)))
title('Incident wave')
xlabel('X');   ylabel('Y');   zlabel('Z');
hold off
view(0,10)


%%% PREPARE OPERATORS
disp('~~~~~~~~~~~~~ PREPARE OPERATORS ~~~~~~~~~~~~~')

% Green kernel function --> G(x,y) = exp(ik|x-y|)/|x-y| 
Gxy         = @(X,Y) femGreenKernel(X,Y,'[exp(ikr)/r]',k);
gradyGxy{1} = @(X,Y) femGreenKernel(X,Y,'grady[exp(ikr)/r]1',k);
gradyGxy{2} = @(X,Y) femGreenKernel(X,Y,'grady[exp(ikr)/r]2',k);
gradyGxy{3} = @(X,Y) femGreenKernel(X,Y,'grady[exp(ikr)/r]3',k);

% Finite elements
u = fem(sphere,typ);
v = fem(sphere,typ);

% Coupling coeff
beta = 1i*k*0.5;

% Finite element mass matrix --> \int_Sx psi(x)' psi(x) dx
Id = integral(sigma,u,v);

% Finite element boundary operator --> \int_Sx \int_Sy psi(x)' G(x,y) psi(y) dx dy 
tic
Sr = 1/(4*pi) .* regularize(sigma,sigma,u,'[1/r]',v);
S  = 1/(4*pi) .* integral(sigma,sigma,u,Gxy,v,tol) + Sr;
toc

% Finite element boundary operator --> \int_Sx \int_Sy psi(x)' dny G(x,y) psi(y) dx dy 
tic
Dr = 1/(4*pi) .* regularize(sigma,sigma,u,'grady[1/r]',ntimes(v));
D  = 1/(4*pi) .* integral(sigma,sigma,u,gradyGxy,ntimes(v),tol) + Dr;
toc

% Final operator : [1i*k*beta*S - (Id/2 + D)]
tic
LHS  = beta.*S - (0.5*Id + D);
toc

% Finite element incident wave trace --> \int_Sx psi(x)' pw(x) dx
RHS = - integral(sigma,u,IW);


%%% SOLVE LINEAR PROBLEM
disp('~~~~~~~~~~~~~ SOLVE LINEAR PROBLEM ~~~~~~~~~~~~~')

% Preconditionneur ILU
tic
[L,U] = ilu(beta.*Sr - (0.5*Id + Dr));
toc

% Solve linear system : [1i*k*beta*S - (Id/2 + D)] = P0
tic
mu = gmres(@(V) LHS*V,RHS,[],tol,100,L,U);
toc

% Jump for derivative
lambda = beta * mu;


%%% INFINITE SOLUTION
disp('~~~~~~~~~~~~~ INFINITE RADIATION ~~~~~~~~~~~~~')

% Plane waves direction
theta = 2*pi/1e3 .* (1:1e3)';
nu    = [sin(theta),zeros(size(theta)),cos(theta)];

% Green kernel function
xdoty        = @(X,Y) X(:,1).*Y(:,1) + X(:,2).*Y(:,2) + X(:,3).*Y(:,3); 
Ginf         = @(X,Y) 1/(4*pi) .* exp(-1i*k*xdoty(X,Y));
gradxGinf{1} = @(X,Y) 1/(4*pi) .* (-1i*k*X(:,1)) .* exp(-1i*k*xdoty(X,Y));
gradxGinf{2} = @(X,Y) 1/(4*pi) .* (-1i*k*X(:,2)) .* exp(-1i*k*xdoty(X,Y));
gradxGinf{3} = @(X,Y) 1/(4*pi) .* (-1i*k*X(:,3)) .* exp(-1i*k*xdoty(X,Y));

% Finite element infinite operators
Sinf = integral(nu,sigma,Ginf,v);
Dinf = integral(nu,sigma,gradxGinf,ntimes(v));
% figure
% spy(Sinf)

% Finite element radiation  
sol = Sinf*lambda - Dinf*mu;

% Analytical solution
ref = sphereHelmholtzGal('inf','dir',1,k,nu,type,X0);
N2 = sprintf( '%.2d', norm(ref-sol,2)/norm(ref,2) )
Ninf = sprintf( '%.2d', norm(ref-sol,'inf')/norm(ref,'inf') )

% Graphical representation
figure
plot(theta,log(abs(sol)),'b',theta,log(abs(ref)),'--r')
title("Relative difference : norm 2 = "+N2+" // norm inf = "+Ninf)


%%% DOMAIN SOLUTION
disp('~~~~~~~~~~~~~ RADIATION ~~~~~~~~~~~~~')

% Green kernel function --> G(x,y) = exp(ik|x-y|)/|x-y| 
Gxy         = @(X,Y) femGreenKernel(X,Y,'[exp(ikr)/r]',k);
gradyGxy{1} = @(X,Y) femGreenKernel(X,Y,'grady[exp(ikr)/r]1',k);
gradyGxy{2} = @(X,Y) femGreenKernel(X,Y,'grady[exp(ikr)/r]2',k);
gradyGxy{3} = @(X,Y) femGreenKernel(X,Y,'grady[exp(ikr)/r]3',k);

% Finite element mass matrix --> \int_Sx psi(x)' psi(x) dx
Id = integral(sigma,u,v);

% Finite element boundary operator --> \int_Sx \int_Sy psi(x)' G(x,y) psi(y) dx dy 
tic
Sbnd = 1/(4*pi) .* (integral(sigma,sigma,u,Gxy,v,tol) + ...
    regularize(sigma,sigma,u,'[1/r]',v));
toc

% Finite element boundary operator --> \int_Sx \int_Sy psi(x)' dny G(x,y) psi(y) dx dy 
tic
Dbnd = 1/(4*pi) .* (integral(sigma,sigma,u,gradyGxy,ntimes(v),tol) + ...
    regularize(sigma,sigma,u,'grady[1/r]',ntimes(v)));
toc

% Boundary solution
Psca = Id\(Sbnd*lambda - (0.5*Id*mu + Dbnd*mu));
Pinc = IW(u.dof);
Pbnd = Pinc + Psca;

% Finite element radiative operator --> \int_Sy G(x,y) psi(y) dy 
tic
Sdom = 1/(4*pi) .* (integral(square.vtx,sigma,Gxy,v,tol) + ...
    regularize(square.vtx,sigma,'[1/r]',v));
toc

% Finite element radiative operator --> \int_Sx \int_Sy psi(x)' grady(G(x,y)) ny.psi(y) dx dy 
tic
Ddom = 1/(4*pi) .* ( integral(square.vtx,sigma,gradyGxy,ntimes(v),tol) + ...
    regularize(square.vtx,sigma,'grady[1/r]',ntimes(v)) );
toc

% Domain solution
Psca = Sdom*lambda - Ddom*mu;
Pinc = IW(square.vtx);
Pdom = Pinc + Psca;


% Annulation sphere interieure
r             = sqrt(sum(square.vtx.^2,2));
Pdom(r<=1.01) = Pinc(r<=1.01);


%%% ANAYTICAL SOLUTIONS FOR COMPARISONS
% Analytical solution
Pbnd1 = sphereHelmholtz('dom','dir',1,k,1.001*sphere.vtx) + IW(sphere.vtx);
Pdom1 = sphereHelmholtz('dom','dir',1,k,square.vtx) + IW(square.vtx);
tic
Pbnd2 = sphereHelmholtzGal('dom','dir',1,k,1.001*sphere.vtx,type,X0) + IW(sphere.vtx);
Pdom2 = sphereHelmholtzGal('dom','dir',1,k,square.vtx,type,X0) + IW(square.vtx);
toc

% Comparison with BEM result
Dbem_21=(Pbnd-Pbnd2)/norm(Pbnd,'inf');
Dbem_22=(Pdom-Pdom2)/norm(Pdom,'inf');
%Comparison with previous analytical solution
Dana_21=(Pbnd1-Pbnd2)/norm(Pbnd1,'inf');
Dana_22=(Pdom1-Pdom2)/norm(Pdom1,'inf');
% Comparison BEM with previous analytical solution
Dgyp_21=(Pbnd-Pbnd1)/norm(Pbnd,'inf');
Dgyp_22=(Pdom-Pdom1)/norm(Pdom,'inf');

%%% GRAPHICAL REPRESENTATIONS
% Total field real part
figure
subplot(1,2,1)
plot(sphere,real(Pbnd))
axis equal;
hold on
plot(square,real(Pdom))
title('Numerical solution real part')
colorbar
hold off
view(0,10)

subplot(1,2,2)
plot(sphere,real(Pbnd2))
axis equal;
hold on
plot(square,real(Pdom2))
title('Analytical solution real part')
colorbar
hold off
view(0,10)

% Total field imag part
figure
subplot(1,2,1)
plot(sphere,imag(Pbnd))
axis equal;
hold on
plot(square,imag(Pdom))
title('Numerical solution imag part')
colorbar
hold off
view(0,10)

subplot(1,2,2)
plot(sphere,imag(Pbnd2))
axis equal;
hold on
plot(square,imag(Pdom2))
title('Analytical solution imag part')
colorbar
hold off
view(0,10)

% Total Field modulus
figure
subplot(1,3,1)
plot(sphere,abs(Pbnd))
axis equal;
hold on
plot(square,abs(Pdom))
title('Numerical solution modulus')
colorbar
hold off
view(0,10)

subplot(1,3,2)
plot(sphere,abs(Pbnd2))
axis equal;
hold on
plot(square,abs(Pdom2))
title('Analytical solution modulus')
colorbar
hold off
view(0,10)

subplot(1,3,3)
plot(sphere,abs(Dbem_21))
axis equal;
hold on
plot(square,abs(Dbem_22))
title('Difference between MyCode & GypsiLab BEM')
colorbar
hold off
view(0,10)


% Norm 2 Comparison
figure
subplot(1,2,1)
plot(sphere,abs(Dana_21))
axis equal;
hold on
plot(square,abs(Dana_22))
title('Difference between current & previous analytic')
colorbar
hold off
view(0,10)

subplot(1,2,2)
plot(sphere,abs(Dgyp_21))
axis equal;
hold on
plot(square,abs(Dgyp_22))
title('Difference between previous analytic & BEM result')
colorbar
hold off
view(0,10)


disp('~~> Michto gypsilab !')


