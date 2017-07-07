% removeCOMMOtion.m:
% Takes a momentum triple [p_1 p_2 p_3] and returns it in the same order
% but with the center of mass motion removed. That is, we are converting
% the momentum vectors from the lab frame to the center-of-momentum (COM)
% frame (sometimes called the molecular frame).
% 
% Inputs:
% * momentum: 9-element row vector containing the momentum triple in
%             the form [px1 py1 ... pz3].
% * masses: 3-element row vector [m1 m2 m3] containing the mass of each
%           atom in atomic mass units [amu].
% 
% Output:
% * out: 9-element row vector containing the momentum triple in the
%        form [px1 py1 ... pz3] in the COM frame.
%
% Notes: * This should only need to be done for experimentally measured
%          momentum vectors. Simulated momentum vectors (e.g. from
%          simulateMomenta.m) are already in the COM frame so this should
%          do nothing to them.
% TODO: * Vectorize this function so it converts an entire set of
%         momentum triples. That way the same function can be used for 1
%         or 100 triples.

function out = removeCOMMotion(momentum, masses)
      % Split the momentum triple into X, Y, Z system components.
      p_X = momentum(1:3:7);
      p_Y = momentum(2:3:8);
      p_Z = momentum(3:3:9);

      massSum = sum(masses);

      % Eliminate center of mass motion.
      p_X = p_X - sum(p_X) .* (masses)/massSum;
      p_Y = p_Y - sum(p_Y) .* (masses)/massSum;
      p_Z = p_Z - sum(p_Z) .* (masses)/massSum;

      % Put the vectors back into the original [p1 p2 p3] form.
      momentum(1:9) = reshape([p_X; p_Y; p_Z], 1, 9);
      
      out = momentum;
end