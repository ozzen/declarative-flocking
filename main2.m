%% Basic Plot
f = figure;

%%
iter = simulation_starting_state;
vel = [vx(iter,:) ;vy(iter,:)];
pos = [x(iter,:) ;y(iter,:)];

minDist = MinPairwiseDistance( x(iter,:), y(iter,:) );

% title(['a:' num2string(params.amax) ', v:' num2string(params.vmax) ', beta:' num2string(params.beta) ', rs:' num2string(params.rs)]);
axis equal
grid on;
ax = gca;
ax.LineWidth = 0.5;
ax.GridColor = [0.7 0.7 0.7];
ax.GridAlpha = 1;
hold on
set(gcf, 'Position', get(0, 'Screensize'))

%%  initPlot;
redBlue = [1 0 0; 0 0 1];
if simplex_policy
    c = redBlue(policy(1,:), :);
else
    c = repmat([1 0 0], params.n, 1);
end
s = 1200*ones(params.n,1); %Marker Area in points-squared

%% Leader Color and Size
if LeaderFollow || predator
    s(end) = 1200;
    c(end,:) = [0 0 0];
end

%%
p = scatter(x(iter,:), y(iter,:), s, c ,'.');

if Speed == 1
    pos_ = pos + v_scale * vel;
    xD = [pos(1,:) ; pos_(1,:)];
    yD = [pos(2,:) ; pos_(2,:)];
    l = plot(xD, yD, 'b', 'LineWidth', 1.3); %velocity plot
elseif Speed == 2
    pos_ = pos + v_scale * vel;
    xD = [pos(1,speed_agent) ; pos_(1,speed_agent)];
    yD = [pos(2,speed_agent) ; pos_(2,speed_agent)];
    l = plot(xD, yD, 'b', 'LineWidth', 1.3); %velocity plot
end

if showCircles
    for i = 1:params.n
        r(i) = rectangle('position', [pos(:,i)' - params.dmin/2 params.dmin params.dmin], 'curvature', [1 1]);
    end
end
%% Add-ons

%Acceleration
if Acceleration == 1
%     xA = [pos(1,:) ; pos(1,:) + 10*(vx(2,:)-vx(1,:))];
%     yA = [pos(2,:) ; pos(2,:) + 10*(vy(2,:)-vy(1,:))];
%     
    xA = [pos(1,:) ; pos(1,:) + Acc_scale * accx(2,:)];
    yA = [pos(2,:) ; pos(2,:) + Acc_scale * accy(2,:)];
    ac = plot(xA,yA,'g', 'LineWidth', 1); %acc plot
    
elseif Acceleration == 2
    xA = [pos(1,Acc_agent) ; pos(1,Acc_agent) + Acc_scale * accx(2,Acc_agent)];
    yA = [pos(2,Acc_agent) ; pos(2,Acc_agent) + Acc_scale * accy(2,Acc_agent)];
    ac = plot(xA,yA,'g', 'LineWidth', 1); %acc plot
end

%Add Obstacles:
if Obstacles
    plotRects(obst_file);
%     [~, cp] = point2rects(pos(:,10)', Rect);
%     cp_h = plot(cp(1), cp(2), 'g*');
end

% Add Target
if Target
    rectangle('Position', [target-1 2 2], 'Curvature', [1,1], 'FaceColor', [0 0.4 0]);
    hold on
end



if pointO
    name = 'pointObstacles25.txt';
    fileID = fopen(name,'r');
    ret = fscanf(fileID,"%f");
    n = ret(1);
    rdata = reshape( ret(2:end) , [2, n] )';
    rdata = rdata(1:n-1,:);
    plot(rdata(:,1),rdata(:,2),'k*', 'MarkerSize', 6);
end


% rect = rectangle('Position',[pos(:,1)' - params.dmin, 2*[params.dmin params.dmin]] ,'Curvature',...
%      [1 1], 'EdgeColor', [0.5 0.5 0.5], 'LineStyle', '--');

if LeaderFollow
    centre = pos(:,1);
else
    centre = sum(pos,2)/params.n;
end

if showRegion
%     hold on
% %    reg = plot_covex_region(traj.linear_A{1}, traj.linear_b{1}, '0', pos(:,end));
%      plot_feasible(traj.linear_A{1}, traj.linear_b{1}, [0;0], [-params.amax; -params.amax], [params.amax; params.amax], ...
% 	      'linecolor', 'b', ...
% 	      'linestyle', '-', ...
% 	      'linesep', 0.4, ...
% 	      'lineangle', 80, ...
% 	      'plot_vertices', 'ko', ...
%           'hold', 1);
    [l_constraints, f_constraints] = plot_convex(f, agnet_number, pos, vel, params);

end

%% Set Window
BoxLength =BoxL * MaxPairwiseDistance(pos(1,:), pos(2,:));
limits = BoxLength/2*[-1 +1 -1 +1];
if FullBox
    height_axis = max(5, max(max(y)) - min(min(y)));
    axis([min(min(x)) max(max(x)) min(min(y))-2.5 min(min(y)) + 1.10 * height_axis]);
    textposy = min(min(y)) + 1.05 * height_axis;
    textposx = min(min(x)) + 0.5 * ( max(max(x)) - min(min(x)) );
    if showText
    textt = text(textposx, textposy, ['Minimum pairwise distance: ' num2str(minDist)], 'FontSize', 15);
    end
else
    axis(centre'*[1 0;1 0;0 1;0 1]' + limits);
    textpos = (7/8) * mean(abs(limits));
    if showText
%     textt = text(centre(1) -15 , centre(2)+textpos, ['Minimum pairwise distance: ' num2str(minDist)], 'FontSize', 15);
      title([ 'Minimum pairwise distance: ' num2str(minDist) ', ' num2str(iter)], 'FontSize', 15);
    end
end

%%
if LeaderFollow
    vbar = vel(:,1)./norm(vel(:,1));
    Vnorm = [1 1]' - dot([1 1]', vbar )*vbar;
    Vnorm = Vnorm/norm(Vnorm);
    tmp = pos(:,1) - 0.2*vel(:,1)./ norm(vel(:,1));
    leaderBarrier = [tmp - 2*Vnorm tmp + 2*Vnorm];
    lb = plot(leaderBarrier(1,:),leaderBarrier(2,:), '--', 'Color', [0 0.5 0]);
end

grid on

%% Plot proximity net
if PlotEdges || pbf
%     for i = 18
    for i = 1:params.n
%         [B, I] = sort(sum((repmat(pos(:,i), 1, Num) - pos).^2,1));
        if PlotEdges == 2 && all(i ~= PlotEdgeId)
            continue;
        end
        for j = 1:params.n
            if (i == j)
                continue;
            end
            if PlotEdges
                edges(i,j) = plot( [pos(1, i), pos(1, j)] , [pos(2, i), pos(2, j)], 'Color',[0 0 0], 'LineWidth', 1.0, 'LineStyle', '-');
%                 set(edges(i,j), 'Visible', 'off');
                
                condition = norm(pos(:,i) - pos(:,j)) > params.rs;
%                 delta_p_ij =  pos(:,i) - pos(:,j);
%                 delta_v_ij =  vel(:,i) - vel(:,j);
%                 condition = barrier_function(delta_p_ij, delta_v_ij, params) > 0;
%                 
                if condition
                    set(edges(i,j), 'Visible', 'off');
                else
                    set(edges(i,j), 'Visible', 'on');
                end
            end
            
            if pbf
                pbf_edges(i,j) = plot( [pos(1, i), pos(1, j)] , [pos(2, i), pos(2, j)], 'Color',[1 0 0], 'LineWidth', 1.5, 'LineStyle', ':');
                pp = [pos(:,i)' ; pos(:,j)'];
                vv = [vel(:,i)' ; vel(:,j)'];
                res = barrierFunction( pp, vv, params.dmin, params.amax );
                disp(['res(' num2str(1) '): ' num2str(res)]);

                if res == 0 || norm(pos(:,i) - pos(:,j)) < params.dmin
                    set(pbf_edges(i,j), 'Visible', 'on'); 
                else
                    set(pbf_edges(i,j), 'Visible', 'off'); 
                end
            end
        end % j
%         if i < I(2)
%             set(edges(i,I(2)), 'Visible', 'on');
%         end
    end % i 
end % if

%% ShowID
if showID
    idtxt = num2cell(1:params.n);
    for i = 1:params.n
        idtxt{i} = num2str(idtxt{i});
    end
    agentIDs = text(pos(1,:)', pos(2,:)'- 0.5, idtxt);
end

F = getframe(f);
writeVideo(myVideo,F);
    
%% Simulation.
Flag = true; %Size change can only happen once.
for i = simulation_starting_state + 1:params.steps-1
    currkey=get(gcf,'CurrentKey');
    if currkey == 'return'
        break
    end
    if mod(i,skip) ~= 0
        continue;
    end
    
    vel = [vx(i,:) ;vy(i,:)];
    pos = [x(i,:) ;y(i,:)];
    
    if PlotEdges || pbf
%         for m = 18
        for m = 1:params.n
%            [B, I] = sort(sum((repmat(pos(:,m), 1, params.n) - pos).^2,1));
            if PlotEdges == 2 && all(m ~= PlotEdgeId)
                continue;
            end
            for j = 1:params.n
 
                if (m == j)
                    continue;
                end
                
                if PlotEdges
                    set(edges(m,j), 'XData', [pos(1, m), pos(1, j)], 'YData', [pos(2, m), pos(2, j)]);
%                     set(edges(m,j), 'Visible', 'off');
                    condition = norm(pos(:,m) - pos(:,j)) > params.rs;
%                     delta_p_ij =  pos(:,m) - pos(:,j);
%                     delta_v_ij =  vel(:,m) - vel(:,j);
%                     condition = barrier_function(delta_p_ij, delta_v_ij, params) > 0;
                    if condition
                        set(edges(m,j), 'Visible', 'off');
                    else
                        set(edges(m,j), 'Visible', 'on');
                    end
                end
                
                if pbf
                    set(pbf_edges(m,j), 'Xdata', [pos(1, m), pos(1, j)], 'Ydata', [pos(2, m), pos(2, j)]);
                    pp = [pos(:,m)' ; pos(:,j)'];
                    vv = [vel(:,m)' ; vel(:,j)'];
                    res = barrierFunction( pp, vv, params.dmin, params.amax );
                    disp(['res(' num2str(i) '): ' num2str(res)]);
                    if res == 0 || norm(pos(:,m) - pos(:,j)) < params.dmin
                        set(pbf_edges(m,j), 'Visible', 'on'); 
                    else
                        set(pbf_edges(m,j), 'Visible', 'off'); 
                    end
                    
                    if norm(pos(:,m) - pos(:,j)) < params.dmin
                        set(pbf_edges(m,j), 'Color', 'g');
                    else
                        set(pbf_edges(m,j), 'Color', 'r');
                    end
                end
                

            end % j
%             if m < I(2)
%                 set(edges(m,I(2)), 'Visible', 'on');
%             end
        end % m
    end % if
    
    if LeaderFollow
        centre = pos(:,1);
    else
        centre = sum(pos,2)/params.n;
        %centre = [0,0]';
    end
    
    pos_ = pos + v_scale * vel;
    
    if Speed == 1
        pos_ = pos + v_scale * vel;
        xD = [pos(1,:) ; pos_(1,:)];
        yD = [pos(2,:) ; pos_(2,:)];
    elseif Speed == 2
        pos_ = pos + v_scale * vel;
        xD = [pos(1,speed_agent) ; pos_(1,speed_agent)];
        yD = [pos(2,speed_agent) ; pos_(2,speed_agent)];
    end
    
    if Acceleration == 1
        %         xA = [pos(1,:) ; pos(1,:) + 10*(vx(i+1,:)-vx(i,:))];
        %         yA = [pos(2,:) ; pos(2,:) + 10*(vy(i+1,:)-vy(i,:))];
        xA = [pos(1,:) ; pos(1,:) + Acc_scale * accx(i+1,:)];
        yA = [pos(2,:) ; pos(2,:) + Acc_scale * accy(i+1,:)];
    elseif Acceleration == 2
        xA = [pos(1,Acc_agent) ; pos(1,Acc_agent) + Acc_scale * accx(i+1,Acc_agent)];
        yA = [pos(2,Acc_agent) ; pos(2,Acc_agent) + Acc_scale * accy(i+1,Acc_agent)];
    end
    
    if Speed == 2
        for indx = 1:numel(speed_agent)
            set(l(indx), 'Xdata', xD(:,indx), 'Ydata', yD(:,indx));
        end
    end
    
    if Acceleration == 2
        for indx = 1:numel(Acc_agent)
            set(ac(indx), 'Xdata', xA(:,indx), 'Ydata', yA(:,indx));
        end
    end
    
    for j = 1:params.n
        if Speed == 1
            set(l(j), 'Xdata', xD(:,j), 'Ydata', yD(:,j));
        end
        if Acceleration == 1
            set(ac(j), 'Xdata', xA(:,j), 'Ydata', yA(:,j));
        end
        if drawTraj == 1
            plot([x(i-1,j) x(i,j)] ,[y(i-1,j) y(i,j)], '-', 'Color', traj_color, 'LineWidth', 2);
        elseif drawTraj == 2
            x_trail = x(max(i-trail_len,1):i,j);
            y_trail = y(max(i-trail_len,1):i,j);
            trail(j) = plot(x_trail, y_trail, '-', 'Color', traj_color, 'LineWidth', 2);
        end
        if showID
            set(agentIDs(j), 'Position', pos(:,j)' - [0, 0.5]);    
        end
        if showCircles
            set(r(j), 'Position', [pos(:,j)' - params.dmin/2 params.dmin params.dmin]);
        end
    end
    if showRegion
%         plot_covex_region(traj.linear_A{i}, traj.linear_b{i}, reg, pos(:,end));
%         plot_feasible(traj.linear_A{1}, traj.linear_b{1}, [0;0], [-params.amax; -params.amax], [params.amax; params.amax], ...
%               'linecolor', 'b', ...
%               'linestyle', '-', ...
%               'linesep', 0.4, ...
%               'lineangle', 80, ...
%               'plot_vertices', 'ko', ...
%               'hold', 1);
        delete(l_constraints);
        delete(f_constraints);
        [l_constraints, f_constraints] = plot_convex(f, agnet_number, pos, vel, params);
    end
    
    if simplex_policy
        c = redBlue(policy(i,:), :);
        set(p, 'Xdata',pos(1,:),'Ydata',pos(2,:), 'CData', c);
    else
        set(p, 'Xdata',pos(1,:),'Ydata',pos(2,:));
    end
    
    %     set(rect, 'Position', [pos(:,1)' - params.dmin, 2*[params.dmin params.dmin]]);
    if LeaderFollow
        vbar = vel(:,1)./norm(vel(:,1));
        Vnorm = [1 1]' - dot([1 1]', vbar )*vbar;
        Vnorm = Vnorm/norm(Vnorm);
        tmp = pos(:,1) - 0.2*vel(:,1)./ norm(vel(:,1));
        leaderBarrier = [tmp - 2*Vnorm tmp + 2*Vnorm];
        set(lb, 'Xdata',leaderBarrier(1,:),'Ydata',leaderBarrier(2,:));
    end
    
    BoxLength = BoxL * MaxPairwiseDistance(pos(1,:), pos(2,:));
    dif = max(BoxLength, max(max(pos') - min(pos')));
    limits = [-dif/2 dif/2 -dif/2 dif/2];
    if ~FullBox
        if true || (MaxPairwiseDistance(pos(1,:), pos(2,:))/BoxLength <= 0.50 || i > 130) && Flag
            BoxLength = BoxL * MaxPairwiseDistance(pos(1,:), pos(2,:));
            limits = BoxLength/2*[-1 +1 -1 +1];
            %Adding a one second pause
    %         for copyframes = 1:myVideo.FrameRate
    %             writeVideo(myVideo,F); 
    %         end
    %         Flag = false;
        end
        axis(centre'*[1 0;1 0;0 1;0 1]' + limits);
        textpos = (7/8) * mean(abs(limits));
        if showText
            minDist = MinPairwiseDistance( x(i,:), y(i,:) );
%             set(textt, 'Position',centre + [-15 textpos]' ,'String', ['Minimum pairwise distance: ' num2str(minDist)]);
            title([ 'Minimum pairwise distance: ' num2str(minDist) ', ' num2str(i)], 'FontSize', 15);
        end
    else
        if showText
        minDist = MinPairwiseDistance( x(i,:), y(i,:) );
        set(textt, 'String', [ 'Minimum pairwise distance: ' num2str(minDist) ', ' num2str(i)]);
        end
    end
    if Obstacles
%         [~, cp] = point2rects(pos(:,10)', Rect);
%         set( cp_h, 'XData', cp(1), 'YData',  cp(2) ); 
    end
    uistack(l, 'top');
    uistack(p, 'top');
    F = getframe(f);
    writeVideo(myVideo,F);
    if drawTraj ~= 0 
        delete(trail);
    end
end

