function [] = plotRects(filename)
    [Rect] = readRects(filename);
    for i = 1 : size(Rect,2)
        hold on
        PlotRect(Rect(i));
    end
end

