function [trajs] = makeTrajs_predator( params , nReps, stateFile, rects, target)
    trajs = [];
    for rept = 1:nReps
        [x, y, vx, vy, ax, ay] = read_data(rept, params.N, stateFile);
%         [x, y, vx, vy, ax, ay] = read_shouvik(stateFile, rept);
        
        jmp = uint16(params.ct / params.dt);
        x = x(1:jmp:end, :);
        y = y(1:jmp:end, :);
        vx = vx(1:jmp:end, :);
        vy = vy(1:jmp:end, :);
        ax = ax(1:jmp:end, :);
        ay = ay(1:jmp:end, :);
        steps = size(x,1);
%         steps = params.steps;
        traj.q = zeros(steps,params.N,params.nd);
        traj.p = zeros(steps,params.N,params.nd);
        traj.a = zeros(steps,params.N,params.nd);
        traj.obstacle = zeros(steps,params.N,params.nd);
        traj.pnet = zeros(steps,params.N-1,params.N-1);
        traj.connectivity = zeros(steps, 1);
        traj.v_converg = zeros(steps, 1);
        traj.irreg = zeros(steps, 1);
        traj.diameter = zeros(steps, 1);
        traj.min_pw_distance = zeros(steps, 1);
        traj.min_obstacle_distance = zeros(steps, 1);
        traj.target = target;
        traj.min_predator_distance = zeros(steps, 1);
        traj.number_of_pairwise_collisions = zeros(steps, 1);
        traj.number_of_predator_collisions = zeros(steps, 1);
       

        for t = 1:steps
            q = [x(t,:)' y(t,:)'];
            p = [vx(t,:)' vy(t,:)'];
            a = [ax(t,:)' ay(t,:)'];

            traj.q(t,:,:) = q;
            traj.p(t,:,:) = p;
            traj.a(t,:,:) = a;
            
            pnet = proximityNet(q(1:end-1,:), params.rad_sensing);
%             pnet = ones(params.N-1);
            traj.pnet(t,:,:) = pnet;
            
            traj.connectivity(t) = connectedComponents(pnet);
            traj.irreg(t) = newIrregularity(q, pnet);

            traj.v_converg(t) = velocityConvergenceNew(p(1:end-1,:), pnet);
            traj.diameter(t) = componentWiseDiameter(q(1:end-1,:), pnet);
            traj.min_pw_distance(t) = MinPairwiseDistance(x(t,1:end-1), y(t,1:end-1));
%             [traj.min_obstacle_distance(t), traj.obstacle(t, :, :)] = minFlockObstacleDistance(q, rects);
            traj.number_of_pairwise_collisions(t) = violation_step(q, params);
            [traj.min_predator_distance(t), traj.number_of_predator_collisions(t)] = minFlockPredatorDistance(q);

        end
        trajs = [trajs,traj];
        if ~mod(rept, 5)
            disp(['Simulations Done:' num2str(rept)]);
        end
    end
end

