function [aprimeIndexes,aprimeProbs] = interp_toolkit(aprimeVals,a_grid)

%aprimeVals=reshape(aprimeVals,[1,N_d*N_a]);
aprimeVals=aprimeVals';%,[1,N_d*N_a]);

a_griddiff=a_grid(2:end)-a_grid(1:end-1); % Distance between point and the next point

temp=a_grid-aprimeVals;
temp(temp>0)=1; % Equals 1 when a_grid is greater than aprimeVals

[~,aprimeIndexes]=max(temp,[],1); % Keep the dimension corresponding to aprimeVals, minimize over the a_grid dimension
% Note, this is going to find the 'first' grid point such that aprimeVals is smaller than or equal to that grid point
% This is the 'upper' grid point

% Switch to lower grid point index
aprimeIndexes=aprimeIndexes-1;
aprimeIndexes(aprimeIndexes==0)=1;

% Now, find the probabilities
aprime_residual=aprimeVals'-a_grid(aprimeIndexes);
% Probability of the 'lower' points
aprimeProbs=1-aprime_residual./a_griddiff(aprimeIndexes);

% Those points which tried to leave the top of the grid have probability 1 of the 'upper' point (0 of lower point)
offTopOfGrid=(aprimeVals>=a_grid(end));
aprimeProbs(offTopOfGrid)=0;
% Those points which tried to leave the bottom of the grid have probability 0 of the 'upper' point (1 of lower point)
offBottomOfGrid=(aprimeVals<=a_grid(1));
aprimeProbs(offBottomOfGrid)=1;


end %end function