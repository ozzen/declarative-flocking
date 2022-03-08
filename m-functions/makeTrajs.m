function [ trajs ] = makeTrajs( params , nReps, stateFile, rects, target)
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
        traj.pnet = zeros(steps,params.N,params.N);
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
        traj.number_of_obstacle_collisions = zeros(steps, 1);
        traj.avg_number_neighbors = zeros(steps, 1);

       

        for t = 1:steps
            q = [x(t,:)' y(t,:)'];
            p = [vx(t,:)' vy(t,:)'];
            a = [ax(t,:)' ay(t,:)'];

            traj.q(t,:,:) = q;
            traj.p(t,:,:) = p;
            traj.a(t,:,:) = a;
            
            pnet = proximityNet(q, params.rad_sensing);
%             pnet = ones(params.N, params.N);
            traj.pnet(t,:,:) = pnet;
            
            traj.connectivity(t) = connectedComponents(pnet);
            traj.irreg(t) = newIrregularity(q, pnet);

%             traj.v_converg(t) = velocityConvergenceNew(p, pnet);
            traj.v_converg(t) = velocityConvegenceGlobal(p, params);
%             traj.diameter(t) = componentWiseDiameter(q, pnet);
            traj.diameter(t) = MaxPairwiseDistance( x(t,:), y(t,:) );
            traj.min_pw_distance(t) = MinPairwiseDistance(x(t,:), y(t,:));
            
%             [traj.min_obstacle_distance(t),...
%              traj.obstacle(t, :, :),...
             traj.number_of_obstacle_collisions(t) = minFlockObstacleDistance(q, rects);
         
            traj.number_of_pairwise_collisions(t) = violation_step(q, params);
            traj.avg_number_neighbors(t) = numberOfNeighbors(pnet);

        end
        trajs = [trajs,traj];
        if ~mod(rept, 5)
            disp(['Simulations Done:' num2str(rept)]);
        end
    end
end

