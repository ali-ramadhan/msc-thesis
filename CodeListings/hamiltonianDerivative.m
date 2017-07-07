% hamiltonianDerivative.m:
% Computes the derivatives of the Hamiltonian (Hamilton's equations)
% for a triatomic molecule.
%
% Inputs:
% * time: a 1x2 row vector containing [InitialTime, FinalTime]
% * p   : a 1x18 row vector containing position and momemtum parameters
%         of the three particles, given as 
%         [x1 y1 z1 x2 y2 z2 x3 y3 z3 px1 px2 ... pz3].
%         p stands for parameters.
%
% Output:
% * out : an nx19 matrix where each row contains [Time, Position[1x9],
%         Momentum[1x9]] in the same format as p. In practice, only
%         the final row is utilized to evaluate the final conditions
%         of the system, i.e. we are only interested in the asymptotic
%         momentum vectors.
%
% Notes: * All units are SI.
%        * This function is written in the form of a system of
%          first-order ordinary differential equations (ODE's),
%          y' = f(t, y), so that it may be solved using MATLAB's
%          4/5th-order Runge-Kutta ODE solver, ode45. Pretty much all
%          numerical ODE solvers expect the system of ODE's to be
%          described in this form.

function out = hamiltonianDerivative(time, p)

% Physical constants
amu = 1.66053886e-27; % [kg], 1 atomic mass unit
e   = 1.60217646e-19; % [C], 1 elementary charge
k   = 8.987551e9;     % [N m^2 C^-2], electrostatic constant

% Atomic masses and charges
m1 = amu*p(19); m2 = amu*p(20); m3 = amu*p(21);
q1 = e*p(22);   q2 = e*p(23);   q3 = e*p(24);

% Calculate the distance between ions. Note that this quantity does
% not preserve vector direction.
r12 = ((p(1)-p(4))^2 + (p(2)-p(5))^2 + (p(3)-p(6))^2)^0.5; % [m]
r13 = ((p(1)-p(7))^2 + (p(2)-p(8))^2 + (p(3)-p(9))^2)^0.5; % [m]
r23 = ((p(4)-p(7))^2 + (p(5)-p(8))^2 + (p(6)-p(9))^2)^0.5; % [m]

% pDot is a column vector with components [vx1; vy1; vz1; vx2; ... vz3;
% p'x1; ... p'z3]. These quantities are produced by taking the first
% derivative of the Hamiltonian with respect to the appropriate variable.

pDot = [p(10)./m1; p(11)./m1; p(12)./m1; ...
        p(13)./m2; p(14)./m2; p(15)./m2; ...
        p(16)./m3; p(17)./m3; p(18)./m3; ...

        k*q1*q2*(p(1)-p(4))/r12^3 + k*q1*q3*(p(1)-p(7))/r13^3; ...
        k*q1*q2*(p(2)-p(5))/r12^3 + k*q1*q3*(p(2)-p(8))/r13^3; ...
        k*q1*q2*(p(3)-p(6))/r12^3 + k*q1*q3*(p(3)-p(9))/r13^3; ...

        k*q2*q1*(p(4)-p(1))/r12^3 + k*q2*q3*(p(4)-p(7))/r23^3; ...
        k*q2*q1*(p(5)-p(2))/r12^3 + k*q2*q3*(p(5)-p(8))/r23^3; ...
        k*q2*q1*(p(6)-p(3))/r12^3 + k*q2*q3*(p(6)-p(9))/r23^3; ...

        k*q3*q1*(p(7)-p(1))/r13^3 + k*q3*q2*(p(7)-p(4))/r23^3; ...
        k*q3*q1*(p(8)-p(2))/r13^3 + k*q3*q2*(p(8)-p(5))/r23^3; ...
        k*q3*q1*(p(9)-p(3))/r13^3 + k*q3*q2*(p(9)-p(6))/r23^3];

out = [pDot; p(19:24)];
end