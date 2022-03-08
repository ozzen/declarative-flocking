function [] = setUpdates( ball, p, l, rect, ac, textt, pos, prevPos, vel, acc, limits, i, params, options)
    pos_ = pos + params.dt * vel;
    xD = [pos(1,:) ; pos_(1,:)];
    yD = [pos(2,:) ; pos_(2,:)];
    Ds = params.Ds;

    if options.Acc
        xA = [ pos(1,:); pos(1,:) + 1 * acc(1,:) ];
        yA = [ pos(2,:); pos(2,:) + 1 * acc(2,:) ];
    end
    for j = 1:params.n
%         set(l(j), 'Xdata', xD(:,j), 'Ydata', yD(:,j)); %Plot vel vector
        set(l(j), 'XData', pos(1,j), 'YData', pos(2,j), 'UData', vel(1,j), 'VData', vel(2,j));
        set(ball(j), 'position', [pos(:,j)' - 0.5*[Ds Ds] [Ds Ds]]);
        if options.Acc
%             set(ac(j), 'Xdata', xA(:,j), 'Ydata', yA(:,j));%Plot acc vector
            set(ac(j), 'XData', pos(1,j), 'YData', pos(2,j), 'UData', acc(1,j), 'VData', acc(2,j));
        end
        plot([prevPos(1,j) pos(1,j)], [prevPos(2,j) pos(2,j)], 'Color', [0.6 0.6 0.6], 'LineStyle', '--', 'LineWidth', 1.4);
    end
    set(p, 'Xdata',pos(1,:),'Ydata',pos(2,:));
%     set(rect, 'Position', [pos(:,1)' - [params.Ds params.Ds] 2*params.Ds 2*params.Ds] );
    
    axis(setAxis( pos, params, limits ));
    textpos = (7/8) * mean(abs(limits));
    minDist = MinPairwiseDistance(pos(1,:), pos(2,:) );
    set(textt, 'Position',sum(pos,2)/params.n + [0 textpos]' ,'String', [num2str(i) ' ' num2str(minDist)]);

end

