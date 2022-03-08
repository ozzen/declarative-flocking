%% view_3D.m is called by Swarms New/main_view_animation.m
% Display video for 3D flock for
%   1. 3D point model dynamics
%   2. Quadcopter model dynamics (6-DOF)
%
%
f = figure;
%% init plot

iter = 1;
vel = [vx(iter,:) ;vy(iter,:); vz(iter,:)];
pos = [x(iter,:) ;y(iter,:); z(iter,:)];

pos_ = pos + v_scale * vel;
xD = [pos(1,:) ; pos_(1,:)];
yD = [pos(2,:) ; pos_(2,:)];
zD = [pos(3,:) ; pos_(3,:)];

c = repmat([1 0 0], params.n, 1);
if params.turn
    for i = traj.edge_agents
        c(i,:) = [0,0,0];
    end
end
sze = 300*ones(params.n,1); %Marker Area in points-squared

rot_speed = [0,0,0,0]';
if ~showQuad
    p = scatter3(x(iter,:), y(iter,:), z(iter,:), sze, c ,'.');
    hold on
    l = plot3(xD, yD, zD, 'b', 'LineWidth', 1.3); %velocity plot
else
    for k = 1:params.n
        p = 0;
        l = 0;
        if rotor_speed
            rot_speed = omegas(iter,:,k)' - mean(omegas(iter,:,k));
        end
        [pn{k}, ln{k}] = display_quadcopter(s(1,:,k), rot_speed, params, p, l, 0);
    end
end
%Acceleration
if Acceleration
    xA = [pos(1,:) ; pos(1,:) + Acc_scale * accx(2,:)];
    yA = [pos(2,:) ; pos(2,:) + Acc_scale * accy(2,:)];
    zA = [pos(3,:) ; pos(3,:) + Acc_scale * accz(2,:)];
    
    ac = plot3(xA,yA,zA,'g', 'LineWidth', 1); %acc plot
end
%text with quadrotors
if isText
    str = {'1','2','3','4','5'};
    txt = text(x(iter,:),y(iter,:),z(iter,:)-1,str);
end

camproj('perspective');
%% Set Window

axis equal
grid on;
ax = gca;
ax.LineWidth = 0.5;
ax.GridColor = [0.7 0.7 0.7];
ax.GridAlpha = 1;
hold on
set(gcf, 'Position', get(0, 'Screensize'))
limits = BoxLength/2*[-1 +1 -1 +1 -1 +1];
if FullBox
    min_x = min(min(x));
    max_x = max(max(x));
    
    min_y = min(min(y));
    max_y = max(max(y));
    
    min_z = min(min(z));
    max_z = max(max(z));
    cush = 0.15;
    %     h_x = 1.15*max(5, max(max(x)) - min(min(x)));
    %     h_y = 1.15*max(5, max(max(y)) - min(min(y)));
    %     h_z = 1.15*max(5, max(max(z)) - min(min(z)));
    h_x = max_x - min_x;
    h_y = max_y - min_y;
    h_z = max_z - min_z;
    
%     axis([min(min(x))-cush max(max(x))+h_x min(min(y))-cush min(min(y))+h_y min(min(z))-cush min(min(z))+h_z]);
    axis([min_x-cush min_x+h_x+cush...
          min_y-cush min_y+h_y+cush...
          min_z-cush min_z+h_z+cush]);
    
else
    centre = mean(pos,2);
    axis(centre'*[1 0 0;1 0 0;0 1 0;0 1 0;0 0 1;0 0 1]' + limits);
end

title(['step number:' num2str(iter)]);
view(azi,ele);
F = getframe(f);
writeVideo(myVideo,F);
%% Simulation
for i = iter+1:params.steps-1
    currkey=get(gcf,'CurrentKey');
    if currkey == 'return'
        break
    end
    
    if mod(i,skip) ~= 0
        continue;
    end
    %     view(45*i/params.steps,45*i/params.steps);
    
    %     azi = azi + turns * (360)/params.steps;
    view(azi,ele);
    
    vel = [vx(i,:); vy(i,:); vz(i,:)];
    pos = [x(i,:); y(i,:); z(i,:)];
    pos_ = pos + v_scale * vel;
    xD = [pos(1,:) ; pos_(1,:)];
    yD = [pos(2,:) ; pos_(2,:)];
    zD = [pos(3,:) ; pos_(3,:)];
    
    if ~showQuad
        set(p, 'Xdata',pos(1,:), 'Ydata',pos(2,:), 'Zdata',pos(3,:));
    end
    
    if Acceleration
        xA = [pos(1,:) ; pos(1,:) + Acc_scale * accx(i+1,:)];
        yA = [pos(2,:) ; pos(2,:) + Acc_scale * accy(i+1,:)];
        zA = [pos(3,:) ; pos(3,:) + Acc_scale * accz(i+1,:)];
    end
    
    for j = 1:params.n
        if ~showQuad
            set(l(j), 'Xdata', xD(:,j), 'Ydata', yD(:,j), 'Zdata', zD(:,j));
        else
            if rotor_speed
                rot_speed = omegas(i,:,j)' - mean(omegas(i,:,j));
            end
            display_quadcopter(s(i,:,j), rot_speed, params, pn{j}, ln{j}, 1);
        end
        if drawTraj == 1
            plot3([x(i-1,j) x(i,j)], [y(i-1,j) y(i,j)], [z(i-1,j) z(i,j)], '-', 'Color', [0.8,0.8,0.8], 'LineWidth', 1.3);
        elseif drawTraj == 2
            x_trail = x(max(i-trail_len,1):i,j);
            y_trail = y(max(i-trail_len,1):i,j);
            z_trail = z(max(i-trail_len,1):i,j);
            trail(j) = plot3(x_trail, y_trail, z_trail, '-', 'Color', [0.7,0.7,0.7], 'LineWidth', 1.3);
        end
        if Acceleration
            set(ac(j), 'Xdata', xA(:,j), 'Ydata', yA(:,j), 'Zdata', zA(:,j));
        end
        if isText
            for tt = 1:params.n
                set(txt(tt), 'Position', [x(i,tt), y(i,tt), z(i,tt)-1]);
            end
        end
    end
    if ~FullBox
        centre = mean(pos,2);
        axis(centre'*[1 0 0;1 0 0;0 1 0;0 1 0;0 0 1;0 0 1]' + limits);
    end
    
    title(['step number:' num2str(i)]);
    F = getframe(f);
    writeVideo(myVideo,F);
    delete(trail);
end

%% plot trajectory
for i = 1:params.n
    plot3( s(:,1,i), s(:,2,i), s(:,3,i) , 'LineWidth', 1.5, 'Color', [0.7,0.7,0.7]);
end
close(myVideo);
drawnow;
