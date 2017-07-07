% simulateMomenta.m:
% Given a list of geometries, simulate a Coulomb explosion for each and
% for each geometry, return the asymptotic momentum vectors for the
% atomic fragments.
%
% Inputs:
% * geometries: nx3 matrix where each row is of the form
%               [r_12 r_23 theta]. r_12 and r_23 should be given in SI
%               units [m] and theta in [deg].
% * masses:     row vector [m1 m2 m3] with the atomic masses in amu.
% * charges:    row vector [q1 q2 q3] with the atomic charges in units
%               of the elementary charge e. So they should be integers.
%
% Output:
% * out: nx12 matrix where each row contains the molecular parameters
%        r_12, r_23, and theta as well as the asymptotic momentum vectors
%        of the three atomic fragments after the Coulomb explosion.
%        [r_12 r_23 theta p_1 p_2 p_3]
%
% Notes: * All units are SI.
%        * parfor loop requires the Parallel Processing Toolbox.

function out = simulateMomenta(geometries, masses, charges, debug)
    nGeometries = size(geometries, 1);

    if debug
    fprintf('# Simulating asymptotic momenta for %d geometries...\n', ...
        nGeometries);
    fprintf('#     Masses = (%.2f, %.2f, %.2f) [amu]\n', ...
        masses(1), masses(2), masses(3));
    fprintf('#     Charges = (%d, %d, %d) [e]\n', ...
        charges(1), charges(2), charges(3));
    end

    out = zeros(nGeometries, 12);

    r_12  = geometries(:, 1);
    r_23  = geometries(:, 2);
    theta = geometries(:, 3);

    % Place each first atom to the left of central atom, along the
    % y-axis.
    x_1 = -r_12;
    y_1 = zeros(nGeometries, 1);
    z_1 = zeros(nGeometries, 1);

    % Place each central atom at the origin.
    x_2 = zeros(nGeometries, 1);
    y_2 = zeros(nGeometries, 1);
    z_2 = zeros(nGeometries, 1);

    % Place each third atom to the right of the central in the +x/+y
    % quadrant, forming an angle (180-theta) with the x-axis.
    x_3 = r_23 .* cosd(180 - theta);
    y_3 = r_23 .* sind(180 - theta);
    z_3 = zeros(nGeometries, 1);

    % For each geometry, calculate the asymptotic momentum vectors it
    % would produce in a Coulomb explosion.
    parfor i = 1:nGeometries
        g = [x_1(i) y_1(i) z_1(i) x_2(i) y_2(i) z_2(i) x_3(i) y_3(i) ...
             z_3(i)]; % initial positions
        p_0 = zeros(1, 9); % initial momentum
        p = coulombExplode([g p_0], masses, charges);
        out(i,:) = [r_12(i) r_23(i) theta(i) p(1:3) p(4:6) p(7:9)];
        
        % Write progress report to console every 100 simulations.
        % Note: Simulations may not be done in order since loop
        % iterations are executed in parallel in a nondeterministic order
        % (parfor loop).
        if debug && (rem(i, 100) == 0)
            fprintf('Simulated geometry #%d/%d.\n', i, nGeometries);
            drawnow('update');
        end
    end

    % Extract each atom's momentum into 2D vectors in preparation to
    % rotate into our convention.
    p_1 = out(:, 4:5);
    p_2 = out(:, 7:8);
    p_3 = out(:, 10:11);

    % Put each momentum into column vector form so we can use matrix
    % multiplication.
    p_1 = p_1'; p_2 = p_2'; p_3 = p_3';

    % Calculate the angle between the central atom and the +x-axis then
    % rotate the three momentum vectors back towards the origin by that
    % much so that the central atom's momentum vector points along the
    % +x-axis.
    for i = 1:nGeometries
        theta_2x = atan2(p_2(2, i), p_2(1, i));
        
        % Rotation matrix for rotating points in the xy-plane
        % counterclockwise through an angle -theta_2x about the origin.
        R = [cos(-theta_2x) -sin(-theta_2x); ...
             sin(-theta_2x) cos(-theta_2x);];
        
        % Rotate each vector. 
        p_1(:, i) = R*p_1(:, i);
        p_2(:, i) = R*p_2(:, i);
        p_3(:, i) = R*p_3(:, i);
    end

    % Put everything back into row vector form.
    p_1 = p_1'; p_2 = p_2'; p_3 = p_3';

    % Set the z components (and also y in case of carbon) to 0 so they
    % all have exactly the exact same numerical value rather than 0.0000
    % and -0.0000, etc. This makes data analysis and filtering more
    % convenient as we can just check for equality with 0.
    % Note: In our momentum vector convention for triatomic molecules,
    % these components should all be zero anyways.
    p_1(:, 3) = 0;
    p_2(:, 2) = 0; p_2(:, 3) = 0;
    p_3(:, 3) = 0;

    out = [r_12 r_23 theta p_1 p_2 p_3];
end

% coulombExplode:
% Given the initial positions and momentum vectors of three ions forming
% a triatomic molecule, simulate their Coulomb explosion and return the
% asymptotic momentum vectors of the three ions.
%
% Inputs:
%  * initialConditions: 18-element row vector containing the initial
%        positions and momentum of the atoms, in the form
%        [x1_0 y1_0 ... z3_0 px1_0 py1_0 ... pz3_0].
% * masses:  row vector [m1 m2 m3] with the atomic masses in amu.
% * charges: row vector [q1 q2 q3] with the atomic charges in units
%            of the elementary charge e. So they should be integers.
%
% Output:
%  * out: A 9-element row vector containing the asymptotic momentum
%         vector components for each ion in the form
%         [px1 py1 pz1 ... pz3]

function out = coulombExplode(initialConditions, masses, charges)
    % Set some error tolerances and initial step sizes for the ODE solver
    % to use. See: https://www.mathworks.com/help/matlab/ref/odeset.html
    options = odeset('AbsTol', 1e-27, ...   % absolute error tolerance
                     'RelTol', 1e-6, ...    % relative //
                     'InitialStep', 1e-18); % initial step size

    % Solve the first-order ODE system described in
    % hamiltonianDerivative.m from t=0 to t=1e-11, by which time the
    % atoms have long attained their asymptotic values. od45 uses
    % adaptive time steps so the intermediate times won't be uniformly
    % distributed.
    [t,y] = ode45('hamiltonianDerivative', [0 1e-11], ...
                  [initialConditions masses charges], options);

    % We have all the intermediate positions and momentum values but we
    % only care about the asymptotic momentum components, so we return
    % them and ignore everything else.
    out = y(size(t, 1), 10:18);
end
