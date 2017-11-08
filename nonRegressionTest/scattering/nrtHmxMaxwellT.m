%+========================================================================+
%|                                                                        |
%|            This script uses the GYPSILAB toolbox for Matlab            |
%|                                                                        |
%| COPYRIGHT : Matthieu Aussal & Francois Alouges (c) 2015-2017.          |
%| PROPERTY  : Centre de Mathematiques Appliquees, Ecole polytechnique,   |
%| route de Saclay, 91128 Palaiseau, France. All rights reserved.         |
%| LICENCE   : This program is free software, distributed in the hope that|
%| it will be useful, but WITHOUT ANY WARRANTY. Natively, you can use,    |
%| redistribute and/or modify it under the terms of the GNU General Public|
%| License, as published by the Free Software Foundation (version 3 or    |
%| later,  http://www.gnu.org/licenses). For private use, dual licencing  |
%| is available, please contact us to activate a "pay for remove" option. |
%| CONTACT   : matthieu.aussal@polytechnique.edu                          |
%|             francois.alouges@polytechnique.edu                         |
%| WEBSITE   : www.cmap.polytechnique.fr/~aussal/gypsilab                 |
%|                                                                        |
%| Please acknowledge the gypsilab toolbox in programs or publications in |
%| which you use it.                                                      |
%|________________________________________________________________________|
%|   '&`   |                                                              |
%|    #    |   FILE       : nrtHmxMaxwellT.m                              |
%|    #    |   VERSION    : 0.30                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal & Francois Alouges            |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 31.10.2017                                    |
%| ( === ) |   SYNOPSIS   : Solve PEC scatering problem with EFIE         |
%|  `---'  |                                                              |
%+========================================================================+
%
% Domaine Omega de frontiere Sigma et de normale sortante n = n^-
%
% Equation :
% Nabla x E - 1i k H = 0 dans Omega^+
% Nabla x H + 1i k E = 0 dand Omega^+
%
% Formulation integrales physiques :
% J = n x H^+       M = E^+ x n
% E(x) = Einc(x) + [T] J(x) + [K] M(x)
% H(x) = Hinc(x) + [T] J(x) - [K] M(x)
%
% Operateurs integraux :
% [T]i,j = 1i k (\int_Ti \int_Tj G(x,y) Jj(y) dot Ji(x)) +
%     - 1i / k (\int_Ti \int_Tj G(x,y) Nabla_y dot Jj(y) Nabla_x dot Ji(x)
%     dx dy
%
% [nxK]i,j = \int_Tj Ji(x) dot n(x) x (\int_Ti \Nabla_y(G(x,y)) x Jj(y)dy) dx
%
% Noyau de Green : G(r) = exp(1i k |r|)/(4 pi |r|) et |r| = |x-y|
%
% Condition aux limites : conducteur parfait sur la frontiere
% E^+ x n = 0 sur \Sigma
%
% EFIE : [T] J = - Einc_tan sur \Sigma
% MFIE : [Id/2 + nxK] J = nxHinc
% CFIE1 : a EFIE + (1-a) MFIE
% CFIE2 : a EFIE + MFIE

% Cleaning
clear all
close all
clc

% Library path
addpath('../../openDom')
addpath('../../openFem')
addpath('../../openHmx')
addpath('../../openMsh')

% Mise en route du calcul parallele 
% matlabpool; 
% parpool

% Parameters
N   = 1e3
tol = 1e-3
typ = 'RWG'
gss = 3

% Spherical mesh
sphere = mshSphere(N,1);
sigma  = dom(sphere,gss);
figure
plot(sphere)
axis equal

% Frequency adjusted to maximum esge size
stp = sphere.stp;
k   = 1/stp(2);
c   = 299792458;
f   = (k*c)/(2*pi);
disp(['Frequency : ',num2str(f/1e6),' MHz']);

% Incident direction and field
X0 = [0 0 -1]; 
E  = [0 1  0]; % Polarization (+x for Theta-Theta and +y for Phi-Phi)
H  = cross(X0,E);

% Incident Plane wave (electromagnetic field)
PWE{1} = @(X) exp(1i*k*X*X0') * E(1);
PWE{2} = @(X) exp(1i*k*X*X0') * E(2);
PWE{3} = @(X) exp(1i*k*X*X0') * E(3);

% Incident wave representation
plot(sphere,real(PWE{2}(sphere.vtx)))
title('Incident wave')
xlabel('X');   ylabel('Y');   zlabel('Z');
hold off
view(0,10)



%%% SOLVE LINEAR PROBLEM
disp('~~~~~~~~~~~~~ SOLVE LINEAR PROBLEM ~~~~~~~~~~~~~')

% Green kernel function --> G(x,y) = exp(ik|x-y|)/|x-y| 
Gxy = @(X,Y) femGreenKernel(X,Y,'[exp(ikr)/r]',k);

% Finite elements
u = fem(sphere,'RWG');
v = fem(sphere,'RWG');

% Finite element boundary operator
tic
LHS = 1i*k/(4*pi)*integral(sigma, sigma, v, Gxy, u, tol) ...
    - 1i/(4*pi*k)*integral(sigma, sigma, div(v), Gxy, div(u), tol) ;
toc

figure
spy(LHS)

% Regularization
tic
LHS = LHS + 1i*k/(4*pi)*regularize(sigma, sigma, v, '[1/r]', u) ...
      - 1i/(4*pi*k)*regularize(sigma, sigma, div(v), '[1/r]', div(u));
toc

figure
spy(LHS)

% Right hand side
RHS = - integral(sigma,v,PWE);

% Solve linear system 
tic
[Lh,Uh] = lu(LHS);
toc
tic
J = Uh \ (Lh \ RHS);
toc

figure
spy(Lh)



%%% INFINITE SOLUTION
disp('~~~~~~~~~~~~~ INFINITE RADIATION ~~~~~~~~~~~~~')

% Plane waves direction
Ninf  = 1e3;
theta = 2*pi/1e3 .* (1:Ninf)';
nu    = [sin(theta),zeros(size(theta)),cos(theta)];

% Green kernel function
xdoty = @(X,Y) X(:,1).*Y(:,1) + X(:,2).*Y(:,2) + X(:,3).*Y(:,3); 
Ginf  = @(X,Y) exp(-1i*k*xdoty(X,Y));

% Finite element infinite operator --> \int_Sy exp(ik*nu.y) * psi(y) dx
Tinf = integral(nu,sigma,Ginf,v, tol);
sol  = 1i*k/(4*pi)*cross(nu, cross([Tinf{1}*J, Tinf{2}*J, Tinf{3}*J], nu));

% Radiation infinie de reference, convention e^(+ikr)/r
nMax = 100; refInf = zeros(Ninf,1);
if E(1) == 1
    for jj = 1:Ninf
        refInf(jj,:) = sphereMaxwell(1, -f, theta(jj), 0.0, nMax);
    end
else
    for jj = 1:Ninf
        [~,refInf(jj)] = sphereMaxwell(1, -f, theta(jj), pi/2, nMax);
    end
end
refInf = refInf ./ sqrt(4*pi);

% Radiations infinies en X
if E(1) == 1
    sol = sin(theta)'.*sol(:,3) - cos(theta)'.*sol(:,1);
else
    sol = sol(:,2);
end

% Erreur
eL2   = norm(refInf-sol,2)/norm(refInf,2)
eLINF = norm(refInf-sol,'inf')/norm(refInf,'inf')
             
% Representation graphique
figure
subplot(1,2,1)
plot(theta,20*log10(abs(sol)),'b',theta,20*log10(abs(refInf)),'r--')

subplot(1,2,2)
plot(theta,real(sol),'--b', theta,imag(sol),'--r', theta, real(refInf),':b', theta,imag(refInf),':r');
drawnow



disp('~~> Michto gypsilab !')