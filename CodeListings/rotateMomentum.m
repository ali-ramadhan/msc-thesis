% rotateMomentum.m:
% Rotate the asymptotic momentum vectors produced by a Coulomb explosion
% such they lie in the xy-plane with the middle atom's momentum vector
% pointing along the +x-axis.
% 
% Input:
% * momenta: 1x9 row vector of momentum components in the form
%            [p2x p2y p2z p1x p1y p1z p3x p3y p3z]
%            WARNING: MIDDLE ATOM FIRST, THEN FIRST, THEN THIRD ATOM!
%
% Output:
% * out: 1x11 row vector containing the momentum triples as well as the
%        theta_v and chi angles in the form
%        [p2x p2y p2z p1x p1y p1z p3x p3y p3z theta_v chi]
%
% Notes: * THIS FUNCTION REQUIRES THE MIDDLE ATOM'S MOMENTUM FIRST, THEN
%          THE FIRST, THEN THIRD! I have no idea why I haven't changed
%          it yet...
%
% TODO: * Vectorize this function so it rotates an entire set of momentum
%         triples. That way the same function can be used for 1 or 100
%         triples.

function out = rotateMomentum (momenta)
    % We want to put the momentum components in the form
    % [px1 py1 pz1 px2 py2 pz2 px3 py3 pz3].
    momenta = reshape(momenta, 3, 3)';

    % If all the vectors are already in the xy-plane, just return them.
    z_components = momenta(:,3);
    if z_components == [0; 0; 0]
        % MATLAB Editor may complain about this equality above but it
        % works for momentum vectors produced by simulateMomenta.m or
        % vectors already rotated using this function.
        % TODO: Assume the zeros may be floating-point and very close to
        % zero.
        
        theta_v = acos(dot(momenta(2,:), momenta(3,:)) / ...
                       norm(momenta(2,:)) / norm(momenta(3,:)));
        chi = acos(dot(momenta(1,:), momenta(2,:)-momenta(3,:))/...
                   norm(momenta(1,:))/norm(momenta(2,:)-momenta(3,:)));
        out =  [momenta(2,:), momenta(1,:), momenta(3,:), theta_v, chi];
        return
    end

    % First, check that all three momentum vectors form a plane.
    % If this is the case, the determinant of the momentum vectors will
    % be zero (or very close).
    
    if det(momenta) < 1e-50
        % The normal vector of the plane is the cross product of two of
        % the vectors.
        normal = cross(momenta(1, :), momenta(2, :));
        
        % Normalise the normal vector!
        normal = normal / norm(normal);
        
        % The components of the normal vector (A,B,C) define the plane.
        % By construction, the plane goes through the origin.  We must
        % now find the angle this plane makes with the xy-plane, which is
        % defined by the normal vector normal_xy = [0 0 1]. This angle is
        % called the dihedral angle.
        normal_xy = [0 0 1];
        dihedral = acos(dot(normal, normal_xy)); % [rad]
        
        % The two planes intersect in a line defined by crossing the two
        % normal vectors.
        intersection = cross(normal, normal_xy);
        
        % if intersection - zeros(1,3) ~= 0 % this is only the case if
        intersection = intersection / norm(intersection);
        
        % Insert all of the relevant information into the matrix which
        % performs the plane rotation.
        momenta = rotatePlane(momenta, dihedral, intersection);
        
        theta_v = acos(dot(momenta(2,:), momenta(3,:))/...
                       norm(momenta(2,:))/norm(momenta(3,:)));
        
        % Determine phi, the minimum angle between one of the terminal
        % atoms and the x-axis.
        % phi = abs(atan2(momenta(2,2), momenta(2,1)));
        
        % For some reason I couldn't get atan2 to work, so let's go with
        % this.
        phi = abs(atan(momenta(2,2)/momenta(2,1)));
        
        if momenta(2,1) >= 0 && momenta(2,2) >= 0    % first quadrant
            phi = phi;
        elseif momenta(2,1) < 0 && momenta(2,2) >= 0 % second quadrant
            phi = pi - phi;
        elseif momenta(2,1) < 0 && momenta(2,2) < 0  % third quadrant
            phi = pi + phi;
        else % fourth quadrant
            phi = 2*pi - phi;
        end

        % Rotate everything (clockwise) through the angle phi
        M = [cos(-phi), -sin(-phi), 0; ...
             sin(-phi), cos(-phi),  0; ...
             0, 0, 1];
        
        momenta(1,:) = (M*(momenta(1,:)'))';
        momenta(2,:) = (M*(momenta(2,:)'))';
        momenta(3,:) = (M*(momenta(3,:)'))';

        % Flip in the y-axis (if necessary) such that the second end
        % atom sits in the +x half plane.
        if momenta(3,2) < 0
            momenta(1:2:3,2) = -momenta(1:2:3,2);
        end
        
        if rand(1) > 0.5
            chi = acos(dot(momenta(1,:), momenta(2,:)-momenta(3,:))/...
                  norm(momenta(1,:))/norm(momenta(2,:)-momenta(3,:)));
        else
            chi = acos(dot(momenta(1,:), momenta(3,:)-momenta(2,:))/...
                  norm(momenta(1,:))/norm(momenta(3,:)-momenta(2,:)));
        end
        
        out = [momenta(2,:), momenta(1,:), momenta(3,:), theta_v, chi];
        
    else
        error('Momentum vectors cannot form a plane!')
    end
end

% rotatePlane:
% Rotate the momentum vectors which define the plane Ax + By + Cz = 0 by 
% the angle dihedral in the line of intersection defined by the direction
% cosines a, b, and c made with the plane z=0. This makes the momentum
% vectors exist only in the x-y plane while retaining their configuration
% and magnitude. The dihedral angle is in radians.
function out = rotatePlane (momenta, dihedral, intersection)
  lineMagnitude = norm(intersection);
  a = intersection(1) / lineMagnitude;
  b = intersection(2) / lineMagnitude;
  c = intersection(3) / lineMagnitude;

  M = [a^2*(1-cos(dihedral))+cos(dihedral), ...
       a*b*(1-cos(dihedral))-c*sin(dihedral), ...
       a*c*(1-cos(dihedral))+b*sin(dihedral); ...
       a*b*(1-cos(dihedral))+c*sin(dihedral), ...
       b^2*(1-cos(dihedral))+cos(dihedral), ...
       b*c*(1-cos(dihedral))-a*sin(dihedral); ...
       a*c*(1-cos(dihedral))-b*sin(dihedral), ...
       b*c*(1-cos(dihedral))+a*sin(dihedral), ...
       c^2*(1-cos(dihedral))+cos(dihedral)];

  for i = 1:3
    newVector(i,:) = (M*(momenta(i,:)'))';
  end

  out = newVector;
end
