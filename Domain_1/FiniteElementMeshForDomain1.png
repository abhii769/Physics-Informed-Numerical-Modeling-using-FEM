%% LAB 5 – Potential Flow Over a Bluff Body (Domain 1)
% Abhigyan Das
% Uses PDE Toolbox for geometry and mesh generation
% Solves Laplace equation: -Δφ = 0
% Computes velocity (u,v) and pressure p = -|∇φ|²

clear; clc; close all;

%% --- 1. Define geometry (outer rectangle - square hole)
R1 = [3,4,-30,30,30,-30,-15,-15,15,15]';   % Outer rectangle
R2 = [3,4,-10,10,10,-10,-5,-5,5,5]';       % Inner square hole
gd = [R1, R2];
ns = char('R1','R2')';                     % Transpose is important here!
sf = 'R1 - R2';
[dl, bt] = decsg(gd, sf, ns);


figure;
pdegplot(dl, 'EdgeLabels','on','FaceLabels','on');
axis equal; title('Geometry of Domain 1 (Rectangular Domain with Square Hole)');
xlabel('x'); ylabel('y');

%% --- 2. Mesh generation
[p, e, t] = initmesh(dl, 'Hmax', 2.0);
figure; pdemesh(p, e, t);
axis equal; title('Finite Element Mesh for Domain 1');

%% --- 3. Assemble stiffness matrix and load vector
A = StiffnessAssembler(p, t, 0);
r = zeros(size(p,2),1);

%% --- 4. Apply boundary conditions
phi = zeros(size(p,2),1);
bc_nodes = [];

% Inlet (x = -30): Neumann BC n·∇φ = 1
for k = 1:size(e,2)
    edgeNodes = e(1:2,k);
    xmid = mean(p(1,edgeNodes));
    if abs(xmid + 30) < 1e-3
        r_local = RobinLoadVector2D(p, e(:,k), @(x,y) 1);
        r = r + r_local;
    % Outlet (x = 30): Dirichlet φ = 0
    elseif abs(xmid - 30) < 1e-3
        bc_nodes = [bc_nodes, edgeNodes];
    end
end

bc_nodes = unique(bc_nodes);
free_nodes = setdiff(1:size(p,2), bc_nodes);

%% --- 5. Solve Aφ = r
phi(free_nodes) = A(free_nodes, free_nodes) \ r(free_nodes);
phi(bc_nodes) = 0;

%% --- 6. Plot iso-contours of potential φ
figure;
pdeplot(p, e, t, 'XYData', phi, 'Contour', 'on', 'ColorMap', 'jet');
axis equal; colorbar;
title('Iso-Contours of Potential \phi (Domain 1)');
xlabel('x'); ylabel('y');

%% --- 7. Compute velocity (u,v) using Hatgradient
nt = size(t,2);
u = zeros(1,nt); v = zeros(1,nt);
xmid = zeros(1,nt); ymid = zeros(1,nt);

for k = 1:nt
    local2glb = t(1:3,k);
    x = p(1,local2glb);
    y = p(2,local2glb);
    [area, b, c] = Hatgradient(x, y);
    phi_local = phi(local2glb);
    u(k) = -sum(phi_local .* b);   % u = -∂φ/∂x
    v(k) = -sum(phi_local .* c);   % v = -∂φ/∂y
    xmid(k) = mean(x);
    ymid(k) = mean(y);
end

%% --- 8. Plot velocity vector field
figure;
quiver(xmid, ymid, u, v, 'AutoScale', 'on', 'AutoScaleFactor', 1.2);
axis equal; grid on;
title('Velocity Vector Field (Domain 1)');
xlabel('x'); ylabel('y');

%% --- 9. Compute pressure p = −|∇φ|²
p_field = -(u.^2 + v.^2);

%% --- 10. Plot iso-pressure contours
figure;
scatter(xmid, ymid, 25, p_field, 'filled');
axis equal; colorbar;
title('Iso-Pressure Contours (p = −|∇\phi|²)');
xlabel('x'); ylabel('y');
grid on;

disp('✅ Domain 1 complete: φ, velocity vectors, and pressure plotted.');

%% ===============================================================
%% --- Supporting Functions (Your Original Ones Corrected)
%% ===============================================================

function A = StiffnessAssembler(p,t,a)
np = size(p,2);
nt = size(t,2);
A  = sparse(np,np);
for k = 1:nt
    local2glb = t(1:3,k);
    x = p(1,local2glb); 
    y = p(2,local2glb);
    [area,b,c] = Hatgradient(x,y);
    Ak = (b*b' + c*c')*area;
    A(local2glb,local2glb) = A(local2glb,local2glb) + Ak;
end
end

function r = RobinLoadVector2D(p,e,a)
np = size(p,2); ne = size(e,2);
r  = zeros(np,1);
for i = 1:ne
    local2glb = e(1:2,i);
    x = p(1,local2glb); 
    y = p(2,local2glb);
    len = sqrt((x(1)-x(2))^2 + (y(1)-y(2))^2);
    xc = mean(x); yc = mean(y);
    abar = a(xc,yc);
    rE = abar*[1;1]*len/2;
    r(local2glb) = r(local2glb) + rE;
end
end

function [area,b,c] = Hatgradient(x,y)
area = abs(det([1 1 1; x; y]))/2;
b = [(y(2)-y(3)); (y(3)-y(1)); (y(1)-y(2))]/(2*area);
c = [(x(3)-x(2)); (x(1)-x(3)); (x(2)-x(1))]/(2*area);
end
