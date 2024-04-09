function [x,y] = shapiro_carlos(y0,ind,data)
% shapiro integrates the Shapiro model for 1D flow in space
%
% PROTOTYPE
%
% [x,y] = shapiro(y0,ind,data)
%
% INPUT:
%
% y0[8x1] Initial conditions of the problem y = [M20, V0, a0, T0, rho0, p0, F0, s0]
% ind[8x1] Variation of Independent Variables
%
%   dA_dx(1) = Variation of Area with respect to x
%   dQ_dx(2) = Variation of Heat with respect to x
%   dWx_dx(3) = Variation of Work with respect to x
%   dH_dx(4) = Variation of Energy Function H Defined in the Shapiro
%                Book with respect to x
%   dX_dx(5) = Variation of Drag Force with respect to x
%   dw_dx(6) = Variation of Mass Flow Rate with respect to x
%   dW_dx(7) = Variation of Molar Mass with respect to x
%   dk_dx(8) = Variation of Gamma with respect to x
%
% data[11x1] Data Defining the System
%
%   L0(1) = Initial Length [m]
%   L(2) = Total Length [m]
%   dx(3) = x Increment in Each Integration [m]
%   A(4) = Area [m^2]
%   W(5) = Molar Mass [kg/mol]
%   f(6) = Friction Coefficient
%   D(7) = Mean Hydraulic Diameter [m]
%   h(8) = Enthalpy per Unit of Mass
%   w(9) = Mass Flow Rate of Gas Stream
%   k(10) = Ratio of Specific Heats or Gamma
%   cp(11) = Specific Heat at Constant Pressure
%
% OUTPUT:
%
% x[tx1] Array of Length Discretization
% y[tx9] Solution of the Shapiro Model
%
%   M2(:,1) = Squared Mach Number [-]
%   V(:,2) = Flow Velocity [m/s]
%   a(:,3) = Sound Speed [m/s]
%   T(:,4) = Absolute Temperature [K]
%   rho(:,5) = Density [kg/m^3]
%   p(:,6) = Static Pressure [Pa]
%   F(:,7) = Impulse Function [??]
%   s(:,8) = Entropy per Unit Mass[]
%
% CONTRIBUTORS:
%
% Fernando Aranda
%
% VERSIONS
%
% 2023-03-14: First Version
% 2023-03-15: Second Version
% -------------------------------------------------------------------------
%% DATA EXTRACTION
L0 = data(1);
L = data(2);
dx = data(3);
k = data(10);
%% INTEGRATION SPAN
nspan = ceil((L-L0)/dx);
xspan = linspace(L0,L,nspan);
%% INTEGRATION
%Solution of the Global Shapiro Equations
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
[x, y] = ode113(@(x,y)shapiro_system(x,y,ind,data), xspan, y0, options);
%% SHAPIRO SYSTEM
    function dy = shapiro_system(x,y,ind,data)
    % Shapiro_systems obtains the system of dependent variables according
    % to Shapiro Model
    %
    % PROTOTYPE
    %
    % dy = shapiro_system(x,y,ind,data)
    %
    % INPUT:
    %
    % x = Variable for Integration
    % y = Array of Variables for Integration
    % ind[8x1] Variation of Independent Variables
    %
    %   dA_dx(1) = Variation of Area with respect to x
    %   dQ_dx(2) = Variation of Heat with respect to x
    %   dWx_dx(3) = Variation of Work with respect to x
    %   dH_dx(4) = Variation of Energy Function H Defined in the Shapiro
    %                Book with respect to x
    %   dX_dx(5) = Variation of Drag Force with respect to x
    %   dw_dx(6) = Variation of Mass Flow Rate with respect to x
    %   dW_dx(7) = Variation of Molar Mass with respect to x
    %   dk_dx(8) = Variation of Gamma with respect to x
    %
    % data[11x1] Data Defining the System
    %
    %   L0(1) = Initial Length [m]
    %   L(2) = Total Length [m]
    %   dx(3) = x Increment in Each Integration [m]
    %   A(4) = Area [m^2]
    %   W(5) = Molar Mass [kg/mol]
    %   f(6) = Friction Coefficient
    %   D(7) = Mean Hydraulic Diameter [m]
    %   h(8) = Enthalpy per Unit of Mass
    %   w(9) = Mass Flow Rate of Gas Stream
    %   k(10) = Ratio of Specific Heats or Gamma
    %   cp(11) = Specific Heat at Constant Pressure
    %
    % OUTPUT:
    %
    % dy[8x1] Solution of the Shapiro Model
    %
    %   dM2(:,1) = Equation of Squared Mach Number
    %   V(:,2) = Equation of Flow Velocity
    %   a(:,3) = Equation of Sound Speed
    %   T(:,4) = Equation of Absolute Temperature
    %   rho(:,5) = Equation of Density
    %   p(:,6) = Equation of Static Pressure
    %   F(:,7) = Equation of Impulse Function
    %   s(:,8) = Equation of Entropy per Unit Mass
    %
    % CONTRIBUTORS:
    %
    % Everybody
    %
    % VERSIONS
    %
    % 2023-03-14: First Version
    % 2023-03-15: Second Version
    % -------------------------------------------------------------------------
    %% VARIABLES SETTLING
    M2 = y(1,1);
    V = y(2,1);
    a = y(3,1);
    T = y(4,1);
    rho = y(5,1);
    p = y(6,1);
    F = y(7,1);
    %% INDEPENDENT VARIABLES EXTRACTION
    dA_dx = ind(1);
    dQ_dx = ind(2);
    dWx_dx = ind(3);
    dH_dx = ind(4);
    dX_dx = ind(5);
    dw_dx = ind(6);
    dW_dx = ind(7);
    dk_dx = ind(8);
    %% DATA EXTRACTION
    A = data(4);
    W = data(5);
    f = data(6);
    D = data(7);
    h = data(8);
    w = data(9);
    k = data(10);
    cp = data(11);
    %% SYSTEM MATRIX DEFINITION
    % Values from A.Shapiro The Dynamics and Thermodynamics of Compressible Fluid Flow, Ronald Press Co., 1953
    FF = [
        M2*[ -2*(1 + (k-1)/2 * M2) / (1-M2),   (1+k*M2) / (1-M2),   k*M2*(1 + (k-1)/2 * M2) / (1-M2),2*(1 + k*M2) * (1 + (k-1)/2 * M2) / (1-M2),   -(1+k*M2) / (1-M2),   -1]; % M2
        V*[ -1/(1-M2), 1/(1-M2),   k*M2/(2*(1-M2)), (1+k*M2) /(1-M2),-1/(1-M2),0]; %V
        a*[ ((k-1)/2)/(1-M2), (1-k*M2)/(2*(1-M2)), -(k*(k-1)*M2*M2)/(4*(1-M2)),-((k-1)/2)*M2*(1+k*M2)/(1-M2),(k*M2-1)/(2*(1-M2)),1/2]; %c
        T*[ M2*(k-1)/(1-M2), (1-k*M2)/(1-M2), -(k*(k-1)*M2*M2)/(2*(1-M2)),-(k-1)*M2*(1+k*M2)/(1-M2),M2*(k-1)/(1-M2),0]; % T
        rho*[ M2/(1-M2), -1/(1-M2),   -k*M2/(2*(1-M2)), -M2*(1+k)/(1-M2),1/(1-M2),0]; %rho
        p*[ k*M2/(1-M2),    -k*M2/(1-M2),   -k*M2*(1+ (k-1)*M2)/(2*(1-M2)),  -2*k*M2*(1 + (k-1)/2 * M2) / (1-M2), k*M2/(1-M2),    0]; %p
        F*[1/(1+k*M2), 0, -k*M2/(2*(1+k*M2)), 0, 0, 0]; % Impulse
        cp*[0, 1, (k-1)*M2/2, (k-1)*M2, 0, 0]; % cp (note 2)
        ];
    %% CREATION OF INDEPENDENT VARIABLES
    x = [dA_dx/A,  (dQ_dx - dWx_dx + dH_dx) / (cp*T ),...
        4*f/D + dX_dx/( 0.5*k*p*A*M2 ) - 2*h*dw_dx/w,  dw_dx/w,  dW_dx/W,   dk_dx/k];
    %% SYSTEM OF EQUATIONS
    dy = sum(FF.*x,2);
    end
end