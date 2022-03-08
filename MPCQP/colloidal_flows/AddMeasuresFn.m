function [ trajs_array ] = AddMeasuresFn( trajs )
% Fields added:
%               minPw - Minimum pairwise distance Steps x 1 vector
%               maxPw - Maximum pairwise distance Steps x 1 vector
%
%
    runs = numel(trajs);
    for i = 1:runs
    %     traj = trajMI(:,i);
        n = trajs(i).params.n;
        steps = size(trajs(i).x, 1);
        for t = 1:steps
            trajs(i).minPW(t) = MinPairwiseDistance(trajs(i).x(t,:), trajs(i).y(t,:));
            trajs(i).maxPW(t) = MaxPairwiseDistance(trajs(i).x(t,:), trajs(i).y(t,:));
            
%             q = squeeze(trajs.q(t,:,:));
%             p = squeeze(trajs.p(t,:,:));
            
            q = [trajs(i).x(t,:)' trajs(i).y(t,:)'];
            p = [trajs(i).vx(t,:)' trajs(i).vy(t,:)'];


            pnet = proximityNet(q, trajs(i).params.Ds);
            trajs(i).pnet(t,:,:) = pnet;
            trajs(i).connectivity(t) = connectedComponents(pnet);
            trajs(i).v_converg(t) = velocityConvergence(p);

    %         traj.irreg(t) = paramFreeIrregularity(q, pnet);
            trajs(i).irreg(t) = newIrregularity(q, pnet, trajs.params.Ds);
            trajs(i).diameter(t) = componentWiseDiameter(q, pnet);
        end
        trajs(i).minPW = trajs(i).minPW';
        trajs(i).maxPW = trajs(i).maxPW';
        trajs(i).connectivity = trajs(i).connectivity';
        trajs(i).v_converg = trajs(i).v_converg';
        trajs(i).irreg = trajs(i).irreg';
        trajs(i).diameter = trajs(i).diameter'; 

    end
    trajs(1).params.steps = steps;
    trajs_array = squeeze(struct2cell(trajs));
end

