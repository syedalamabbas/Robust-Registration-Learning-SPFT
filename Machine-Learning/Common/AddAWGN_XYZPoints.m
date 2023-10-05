function [ pointCloudOutput ] = AddAWGN_XYZPoints( pointCloudInput, SNRdB )
%ADDAWGN_XYZPOINTS is the function
pointCloudOutput = zeros(size(pointCloudInput));

isPlotting = 0;

% X - axis
verticesX = pointCloudInput(:,1);
newVerticesX = awgn(verticesX, SNRdB, 'measured');
pointCloudOutput(:,1) = newVerticesX;
if(isPlotting)
    PlotOneDimenNoisySignals(verticesX, newVerticesX);
end


% Y - axis
verticesY = pointCloudInput(:,2);
newVerticesY = awgn(verticesY, SNRdB, 'measured');
pointCloudOutput(:,2) = newVerticesY;
if(isPlotting)
    PlotOneDimenNoisySignals(verticesY, newVerticesY);
end


% Z - axis
verticesZ = pointCloudInput(:,3);
newVerticesZ = awgn(verticesZ, SNRdB, 'measured');
pointCloudOutput(:,3) = newVerticesZ;
if(isPlotting)
    PlotOneDimenNoisySignals(verticesZ, newVerticesZ);
end

    function PlotOneDimenNoisySignals(vertices, newVertices)
        figure, plot(newVertices, 'b')
        hold on
        plot(vertices, 'r')
        hold off
        legend('Noisy points', 'Actual points')
    end
end

