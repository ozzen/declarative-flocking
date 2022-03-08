function [range] = plotLogFileData( fname, yLabel, Num)
A = ReadFile( fname, "%f" );
A = log10(A);
range = [max(A), min(A)]; 
xLabel = 'Time';
plotCurve(A, xLabel, yLabel, 1.5, Num)
grid on
end

