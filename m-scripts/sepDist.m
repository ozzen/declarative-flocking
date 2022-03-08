figure;
for index = 1:2
    xc = x(index,:);
    yc = y(index,:);
    TRI = delaunay(xc,yc);
    TRI = [TRI TRI(:,1)];
    
    Edges = zeros(Num);
    for i = TRI'
        for j = 1:3
            Edges( i(j), i(j+1) ) = 1;
        end
    end
    
    k = Edges == 1;
    % imshow(k '| k )
    k = k'|k;
  
    dist = zeros(3*Num,1);
    add = 0;
    subplot(1,2,1);
    hold on
    for i = 1:Num
        for j = i:Num
            if k(i,j)
                add = add+1;
                dist(add) = norm([xc(i)-xc(j) yc(i)-yc(j)]);
                plot([xc(i) xc(j)], [yc(i) yc(j)], '*--', 'MarkerEdgeColor', 'b', 'Color','r',  'LineWidth', 2 );
            end
        end
    end
    
    axis equal
    dist = dist(1:add);
    subplot(1,2,2);
    histogram(dist, 7);
%     xlim([0 5]);
    hold off
    title(num2str(index));
    pause(0.1);
    clf('reset');
end