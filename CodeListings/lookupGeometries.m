% lookupGeometries.m:
% Takes in a 12-column momentum lookup table `table`, an auxillary lookup
% table `auxtable` for storing results from extra simulations, a list of
% momentum triples `momenta` that you want to find geometries for and an
% error tolerance to aim for.
% 
% It will return a matrix bestGeometries with the best matching
% geometries for each momentum triple in momenta. If one could not be
% found, it will return a row of zeros. It will also return an updated
% auxtable with extra entries computed when running extra simulations.
% 
% Outline of the lookup table implementation:
% 1. Check the momentum triple against all entries in the coarse lookup
%    table and auxillary lookup table, and find the best matching
%    geometry by finding the geometry and momentum triple that minimizes
%    the norm squared of the difference of the two momentum vectors.
% 2. If the geometry we found isn't good enough, simulate the experiment
%    for geometries close to the best one we found from the lookup table
%    and check if the best simulated geometry is good enough. If not,
%    keep doing this iteratively for smaller regions around the best
%    geometry so far until a good enough geometry has been found.
% 3. If successive iterations of step 2 don't yield a good enough
%    geometry, we look at the top 100 geometries we found from step 1
%    and try a different geometry. We try to pick a geometry from a
%    different region in the parameter space in case there is more than
%    one converging region.
% 4. If we have tried converging on three different geometries with no
%    luck, just return zeros. We give up.
% 
% Inputs:
% * table: nx12 matrix with each row containing geometries and their post
%          Coulomb explosion momentum vectors in the form
%          [r12 r23 theta p1x p1y ... p3z]
% * auxtable: optional mx12 matrix. You could put in an empty matrix or
%             an existing momentum lookup table with the same format as
%             `table`. lookupGeometries will use this table in addition
%             to the coarse lookup table but whenever lookupGeometries
%             simulates extra geometries, it will add the simulation to
%             the table so it can be used in future calculations. I.e.
%             memoization, or dynamically updating the table. It's mostly
%             to make sure we don't run the same simulation twice as it's
%             always faster to check against a lookup table than it is to
%             simulate more data.
% * momenta:
% * tolerance: error threshold below which a geometry is considered to be
%              a good, or optimal, solution. Generally, using just the
%              coarse lookup table, you can get down to errors of ~5e-48,
%              which is barely 2 significant figures of accuracy if
%              you're looking at angstroms and degrees. 1e-50 is about
%              2-3 significant figures, and 1e-52 is probably 3-ish. Of
%              course, the lower the desired error tolerance, the longer
%              it takes to converge and sometimes it might not be able
%              to if it gets stuck in a local minimum, in which case a
%              row of zeros is returned.
% 
% Outputs:
% * bestGeometries:
% * auxtable:
% 
% Notes: * All momentum triples should be ordered as OCS in table,
%          auxtable and momenta matrices!
%        * Lookup table can be generated using simulateMomenta.m. You
%          must be careful that the momentum vectors are all in some sort
%          of standard configuration or convention. For example, all OCS
%          momentum vectors have the p_C vector pointing along the
%          +x-axis, the p_S vector always points in the -x,+y direction
%          and the p_O vector points in the -x,-y or +x,-y direction,
%          it's on either side of the -y-axis. All momentum vectors have
%          no z-components in the case of OCS and CO2 as they are
%          triatomic and we can rotate the momentum vectors to lie in a
%          plane. simulateMomenta.m will produce momentum vectors with
%          this configuration and rotateMomentum will put momentum
%          vectors in that configuration if you're processing
%          experimentally measure data.
% 
% TODO: * Have the auxillary table update during simulations.
%       * If it doesn't converge, maybe include the best geometry it
%         found anyways then filter based on fitness?
%       * Multiple debug message levels. Does MATLAB have a logger class?

function [bestGeometries, auxtable] = ...
    lookupGeometries(table, auxtable, momenta, tolerance)

    simulationSteps = 5;
    geometriesToTry = 3;
    coarseR_12Step = 0.05e-10; % inital rStep?
    coarseR_23Step = 0.05e-10; % initial thetaStep?
    coarseThetaStep = 0.25;
    separationThreshold = 10;
    successiveIterationThreshold = 0.05;

    % For each momentum triple you want to find a geometry for...
    for i = 1:size(momenta, 1)
        momentum = momenta(i,:);
        
        % If we're dealing with real data, we want to remove COM (center
        % of mass) motion and rotate the momentum into a standard
        % configuration. If you're dealing with simulated data, that's
        % already been done so you should comment these lines out.
        
        % TODO: only call removeCOMMotion is sum of p is not zero and
        % only call rotateMomentum2 if p vectors aren't in a standard
        % configuration. Easy to check and no need to comment.
        
        %momentum = removeCOMMotion(momentum);
        %momentum = rotateMomentum2(momentum);
        
        % Search through the coarse and auxillary momentum tables for a
        % best matching geometry. Don't expect an error of lower than
        % ~5e-48 from just the coarse lookup table.
        [bestGeometryCoarse, bestErrorCoarse] = ...
            lookupGeometry(table, momentum);
        [bestGeometryAux, bestErrorAux] = ...
            lookupGeometry(auxtable, momentum);
    
        % Save the best geometry we found from the two lookup tables.
        if bestErrorCoarse < bestErrorAux
            bestGeometry = bestGeometryCoarse;
            bestError = bestErrorCoarse;
            fprintf(['[%d] Found coarse geometry (%e, %e, %f) with' ...
                'e = %.2e.\n'], i, bestGeometry(1), bestGeometry(2), ...
                bestGeometry(3), bestErrorCoarse);
        else
            bestGeometry = bestGeometryAux;
            bestError = bestErrorAux;
            fprintf(['[%d] Found aux geometry (%e, %e, %f) with' ...
                'e = %.2e.\n'], i, bestGeometry(1), bestGeometry(2), ...
                bestGeometry(3), bestErrorAux);
        end
        
        % The geometry number we're looking at. It will only go up if we
        % fail to converge on the first coarse geometry we find.
        geometryNum = 1;
        
        topGeometries = top100(table, momentum);
        for g = geometryNum:size(topGeometries,1)
            g2 = topGeometries(g,1:3);
            
            if (g2(1) < 0.95e-10) || (g2(1) > 4.50e-10) || ...
                    (g2(2) < 0.95e-10) || (g2(2) > 4.50e-10)
                continue;
            end
            
            bestGeometry = topGeometries(g,:);
            bestError = norm(bestGeometry(4:12) - momentum)^2;
            fprintf(['[%d] Starting with geometry %d, (%e, %e, %f)' ...
                     'with e = %.2e.\n'], i, g, bestGeometry(1), ...
                     bestGeometry(2), bestGeometry(3), bestError);
                break;
        end
    
        % If we haven't found a good enough geometry, we'll have to run
        % simulations for geometries close to what we converged on to see
        % if we can find a better one. Initially, the lookup table has a
        % resolution of 0.05 A and 0.25 degrees so we start with those
        % then the algorithm will make them smaller as it goes.
        rStep = 0.05e-10;
        thetaStep = 0.25;
        
        originalGeometry = bestGeometry;
        originalError = bestError;
    
        while bestError > tolerance
            rStep = rStep / 5;
            thetaStep = thetaStep / 5;
            
            simulatedTable = simulateRange(bestGeometry, i, rStep, ...
                thetaStep, simulationSteps);
            [bestGeometrySim, bestErrorSim] = ...
                lookupGeometry(simulatedTable, momentum);
            
            if abs(bestErrorSim - bestError)/bestError < 0.05
                % This will happens if we fail to converge onto a
                % geometry so even after simulating a deeper region, we
                % still haven't found anything better. We look at another
                % region if the error is less than 1% better for
                % consecutive geometries. We usually expect the error to
                % be ~1000% better per iteration for a nice convergence
                % point.
            
                if simulationSteps == 5 && geometryNum >= 2
                    fprintf(['[%d] Enhancing momentum lookup table' ...
                        'simulations from %d to %d!\n'], i, 10^3, 20^3);
                    
                    bestGeometry = originalGeometry;
                    bestError = originalError;
                    simulationSteps = 10;
                    rStep = 0.05e-10;
                    thetaStep = 0.25;
                else
                    simulationSteps = 5;
                    
                    if geometryNum == 3
                        fprintf(['[%d] We have looked at %d '...
                            'geometries. Giving up now.\n'], ...
                            i, geometryNum);
                        bestGeometry = zeros(1,12);
                        break;
                    end
                    
                    geometryNum = geometryNum+1;
                    fprintf(['[%d] Cannot converge to a geometry.' ...
                        'Attempting to switch to a better' ...
                        'geometry...\n'], i);
                    
                    topGeometries = top100(table, momentum);
                    
                    space = [topGeometries(:,1:3) ...
                             kron(topGeometries(:,1), false)];
                    g1 = bestGeometrySim(1:3);
                    
                    switched = false;
                    for g = geometryNum:size(topGeometries,1)
                        g2 = topGeometries(g,1:3);
                        
                        if g2(1) < 0.95e-10 || g2(1) > 4.50e-10 || ...
                                g2(2) < 0.95e-10 || g2(2) > 4.50e-10
                            continue;
                        end
                        
                        if ~isConnected(space, [g1 false], [g2 false])...
                                && separation(g1, g2) >= 10
                            bestGeometry = topGeometries(g,:);
                            bestError = ...
                                norm(bestGeometry(4:12) - momentum)^2;
                            fprintf(['[%d] Switched to disconnected' ...
                                'geometry %d, (%e, %e, %f) with' ...
                                'e = %.2e.\n'], i, g, bestGeometry(1),...
                                bestGeometry(2), bestGeometry(3), ...
                                bestErrorSim);
                            switched = true;
                            break;
                        end
                    end
                
                if ~switched
                    for g = geometryNum:size(topGeometries,1)
                        g2 = topGeometries(g,1:3);
                        g2(1:2)
                        
                        if g2(1) < 0.95e-10 || g2(1) > 4.50e-10 || ...
                                g2(2) < 0.95e-10 || g2(2) > 4.50e-10
                            continue;
                        end
                        
                        bestGeometry = topGeometries(g,:);
                        bestError = ...
                            norm(bestGeometry(4:12) - momentum)^2;
                        fprintf(['[%d] Switched to connected geometry'...
                            '%d, (%e, %e, %f) with e = %.2e.\n'], ...
                            i, geometryNum, bestGeometry(1), ...
                            bestGeometry(2), bestGeometry(3), ...
                            bestErrorSim);
                        break;
                    end
                end
                
                originalGeometry = bestGeometry;
                originalError = bestError;
                rStep = 0.05e-10;
                thetaStep = 0.25;
                end
            elseif bestErrorSim < bestError
                bestGeometry = bestGeometrySim;
                bestError = bestErrorSim;
                fprintf(['[%d] Found simulated geometry (%e, %e, %f)' ...
                    'with e = %.2e.\n'], i, bestGeometry(1), ...
                    bestGeometry(2), bestGeometry(3), bestErrorSim);
            end
        end
        simulationSteps = 5;
        
        % This will fail if no simulation was needed.
        fprintf(['[%d] Converged to geometry (%e, %e, %f) with' ...
            'e = %.2e.\n\n'], i, bestGeometry(1), bestGeometry(2), ...
            bestGeometry(3), bestErrorSim);
        
        % Doing 1:12 because sometimes (after enhancing), we end up with
        % a 13th fitness/error column and I'm not sure why.
        bestGeometries(i,:) = bestGeometry(1:12);
    end
end

function [bestGeometry, error_min]  = lookupGeometry(table, momentum)
    error_min = 100;
    j_min = -1;
    
    for j = 1:size(table,1)
        error = norm(momentum - table(j, 4:12))^2;
        
        if error < error_min
            error_min = error;
            j_min = j;
        end
    end
    
    if j_min == -1
        bestGeometry = zeros(12);
    else
        bestGeometry = table(j_min,:);
    end
end

function out = simulateRange(geometry, i, rStep, thetaStep, ...
                             simulationSteps)
    midR12   = geometry(1);
    midR23   = geometry(2);
    midTheta = geometry(3);

    nStep = simulationSteps;
    
    fprintf(['[%d] Simulating data for r_12 = [%e,%e], ' ...
             'r_23 = [%e,%e], theta = [%f,%f].\n'], ...
             i, midR12 - nStep*rStep, midR12 + nStep*rStep, ...
             midR23 - nStep*rStep, midR23 + nStep*rStep, ...
             midTheta - nStep*thetaStep, midTheta + nStep*thetaStep);
    
    j = 1;
    for r_12 = midR12 - nStep*rStep : rStep : midR12 + nStep*rStep
        for r_23 = midR23 - nStep*rStep : rStep : midR23 + nStep*rStep
            for theta = midTheta - nStep*thetaStep : thetaStep : midTheta + nStep*thetaStep
                out(j,:) = simulateMomentum(r_12, r_23, theta);
                j = j+1;
            end
        end
    end
end

function out = top100(table, momentum)
    % Keeping commented because if you use this function, you already
    % have COM removed and momentum vectors rotated.

    %momentum = removeCOMMotion(momentum);
    %momentum = rotateMomentum2(momentum);
    
    out = zeros(size(table,1), 13);
    for i = 1:size(table,1)
        out(i,:) = [table(i,:) norm(momentum - table(i,4:12))^2];
    end

    out = sortrows(out, 13); % Sort rows by error.
    out = out(1:100,:);
end

% space is an nx4 matrix with (r_12, r_23, theta, false) data points
% corresponding to possible geometries from the lookup table. p1 and p2
% are possible geometries. isConnected returns true if p1 and p2 are
% within the same region of convergence, that is, they exists a path
% between p1 and p2 in the lookup table. Otherwise, it returns false.
% isConnected is usually used to check if two points in the top 100
% geometries are in the same region or not. If they are not within the
% same region, then we can use that point to check if it represents the
% geometry we're looking for to avoid searching for the right geometry
% simply locally.
function out = isConnected(space, p1, p2)
    [contained,idx] = ...
        ismember(single(p1(1:3)), single(space(:,1:3)), 'rows');

    % If p1 is not in space, return false. If p1 is in space and we've
    % already checked it, then return false. Only continue if p1 is in
    % space and we haven't already checked it.
    if idx == 0
        % disp('not in table');
        out = false;
        return;
    else
        if space(idx,4)
            % disp('already visited');
            out = false;
            return;
        else
            % disp('new')
            space(idx,4) = true;
        end
    end
    
    if all(single(p1) == single(p2))
        %disp('there!');
        out = true;
        return;
    end

    r_12Step  = 0.05e-10;
    r_23Step  = 0.05e-10;
    thetaStep = 0.25;
    
    P = p2 - p1; % Direction from geometry p1 to p2.
    
    possibleSteps = ...
        [dir(P,1)*r_12Step 0                 0;
         0                 dir(P,2)*r_23Step 0;
         0                 0                 dir(P,3)*thetaStep;
         dir(P,1)*r_12Step dir(P,2)*r_23Step 0;
         dir(P,1)*r_12Step 0                 dir(P,3)*thetaStep;
         0                 dir(P,2)*r_23Step dir(P,3)*thetaStep;
         dir(P,1)*r_12Step dir(P,2)*r_23Step dir(P,3)*thetaStep;];
    
    possibleSteps = unique(possibleSteps, 'rows');
    possibleSteps( ~any(possibleSteps,2), : ) = [];
             
    for i = 1:size(possibleSteps,1)
        if isConnected(space, p1 + [possibleSteps(i,:) false], p2)
            out = true;
            return;
        end
    end
    out = false;
end

function out = dir(v, d)
    if v(d) == 0
        out = 0;
    elseif v(d) > 0
        out = 1;
    else
        out = -1;
    end
end
    
function out = separation(g1, g2)
    r_12Step = 0.05e-10;
    r_23Step = 0.05e-10;
    thetaStep = 0.25;
    
    delta_g = abs(g2 - g1);
    sep = delta_g(1)/r_12Step + delta_g(2)/r_23Step + ...
          delta_g(3)/thetaStep;
    out = sep;
end