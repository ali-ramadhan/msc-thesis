% multiStartTriatomic.m:
% ...
% 
% Inputs:
% * momenta: 
% * masses:
% * charges: 
% * fOutFilenamePrefix: 
% * startingIndex: 
% * runs:
% 
% Output: None.
% 
% Notes: * This function does not return anything as the reconstructions
%          are done in parallel in a parfor loop. So to avoid multiple
%          threads writing to the same out variable, the results of each
%          reconstruction are saved to a separate file called
%          `fOutFilenamePrefix_Gxxxxx.log` where xxxxx is the index of
%          the geometry in the `momenta` matrix. Yes, there will be
%          problems if you're trying to reconstruct more than 99,999
%          geometries...

function multiStartTriatomic(momenta, masses, charges, ...
                             fOutFilenamePrefix, startingIndex, runs)
    nMomenta = size(momenta, 1);

    parfor i = startingIndex:nMomenta
        fOutFilename = strcat(fOutFilenamePrefix, '_G', ...
                              sprintf('%05d',i), '.log');
        fOut = fopen(fOutFilename, 'a');
        
        p = momenta(i, :);
        p = removeCOMMotion(p);
        p = rotateMomentum2(p);
        
        pGoal = p;
        residualNormObjective = @(g)residualNorm(g, pGoal, masses, ...
                                                 charges);

        % These will be our bounds for the multi start algorithm. The
        % algorithm will not search outside of these bounds. They include
        % a wide variety of possible geometries but nothing super
        % unrealistic (e.g. super compressed bonds). Lengths are in [pm]
        % and angles in [deg] because computers don't like numbers that
        % are too small (it seems to be harder to converge when one
        % parameter is ~1e-12 and another ~1e2.) and I thought it would
        % be nice to keep all numbers in the same order of magnitude.
        r_12LowerBound  = 100;
        r_12UpperBound  = 500;
        r_23LowerBound  = 100;
        r_23UpperBound  = 500;
        thetaLowerBound = 140;
        thetaUpperBound = 180;
        lowerBounds = [r_12LowerBound r_23LowerBound thetaLowerBound];
        upperBounds = [r_12UpperBound r_23UpperBound thetaUpperBound];
        
        % You have to give the multi start algorithm a starting point so
        % I thought might as well give it some middle point. It's not
        % sensitive to the starting point at all which is great. It's
        % going to guess different points anyways.
        % Note: r_12 = 115 pm, r_23 = 156 pm, theta = 175 deg is ground
        % state equilibrium.
        r_12Initial  = 250;
        r_23Initial  = 250;
        thetaInitial = 170;
        initialGeometry = [r_12Initial r_23Initial thetaInitial];

        options = optimoptions('fmincon', ...
                               'Algorithm', 'interior-point', ...
                               'Display', 'off', ...
                               'MaxFunEvals', 3000);
                           
        problem = createOptimProblem('fmincon', ...
            'objective', residualNormObjective, ...
            'lb', lowerBounds, 'ub', upperBounds, ...
            'x0', initialGeometry, ...
            'options', options);
        
        ms = MultiStart('UseParallel', 'always', ...
                        'Display', 'off', ...
                        'StartPointsToRun', 'bounds');
        
        [g, fval, exitflag, output, solutions] = run(ms, problem, runs);

        fprintf('G%05d DONE @ %s.\n', i, datestr(now));
        fprintf(['G%05d Best geometry found: ' ...
                 '(%.2f pm, %.2f pm, %.2f deg)' ...
                 'with log residual norm %.2f and exit flag %d.\n'], ...
                 i, g(1), g(2), g(3), fval, exitflag);
        fprintf('G%05d Solver: funcCountNumber:       %d\n', ...
                i, output.funcCount);
        fprintf('G%05d         localSolverIncomplete: %d\n', ...
                i, output.localSolverIncomplete);
        fprintf('G%05d         localSolverNoSolution: %d\n', ...
                i, output.localSolverNoSolution);
        fprintf('G%05d         localSolverSuccess:    %d\n', ...
                i, output.localSolverSuccess);
        fprintf('G%05d         localSolverTotal:      %d\n', ...
                i, output.localSolverTotal);
        
        % We don't always get back 'runs' solutions so we count how many
        % we found.
        numSolutionsFound = size([solutions.Fval], 2);
        
        % We put each distinct solution we found into a row vector and
        % print them all.
        mostLikelyGeometries = ...
            [i*ones(numSolutionsFound, 1), ...
             reshape([solutions.X], 3, numSolutionsFound)', ...
             [solutions.Fval]', ...
             [solutions.Exitflag]'];
        
        fprintf('G%05d Writing %s.\n', i, fOutFilename);
        for j = 1:numSolutionsFound
            fprintf('G%05d G %d\t%3.6f\t%3.6f\t%3.6f\t%2.2f\t%d\n', ...
                    i, mostLikelyGeometries(j,:));
            fprintf(fOut, '%d\t%3.6f\t%3.6f\t%3.6f\t%2.2f\t%d\n', ...
                    mostLikelyGeometries(j,:));
        end
        
        fprintf('\n');
        fprintf(fOut,'\n');
        fclose(fOut);
    end
end

% This is our objective or fitness function for the multi start
% algorithm. It takes a vector g with our molecule's configuration
% (r_12, r_23, theta) and simulates a Coulomb explosion for it. It then
% compares the asymptotic momentum produced with the goal momentum we are
% attempting to recreate and returns the log10 of the norm of the
% difference between the two momentum squared.
function rn = residualNorm(g, pGoal, masses, charges)
    g = [1e-12*g(1) 1e-12*g(2) g(3)];
    p = simulateMomentum(g, masses, charges);
    p = p(4:12);
    rn = log10(norm(pGoal - p)^2);
end

% rotateMomentum2 rotates the momentum to eliminate the z-component. We
% change format 123->213 (e.g. OCS->COS) because rotateMomentum expects
% 213/COS.
function out = rotateMomentum2(momentum)
    p_1 = momentum(1:3);
    p_2 = momentum(4:6);
    p_3 = momentum(7:9);

    momentum = [p_1 p_2 p_3];
    momentum = rotateMomentum(momentum);
    momentum = [momentum(4:6) momentum(1:3) momentum(7:9)];

    out = momentum(1:9);
end
