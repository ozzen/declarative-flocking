function [ ret ] = cost_sum( u, pos, vel, params )
%% cmpc
    N = params.n;
    acc = u2acc(u, params);
    control_steps = uint8(params.ct / params.dt);
    res = 0;
    for h = 1:params.h
        a = acc(:,:,h)';
        for k = 1:control_steps
            if params.predator
                a(:,end) = [0;0]; % predator does not accelerate in the model
%                 [a(:,end)] = controller_predator(pos, params); %predator's accurate model
            end
            [pos, vel] = dynamics(pos, vel, a, params);
        end
        res = res + fitness( pos, vel, params );
    end
    ret = res/params.h + sum(u.^2)/params.h;
end

