
function [ntrpy] = Bentropy (histog)

histog(histog==0) = [];       % eliminating zeros
histog = full(histog);        % sparse to double

probs = histog./sum(histog);

ntrpy = -sum(probs.*log2(probs));

% ntrpy  = 0;
% 
% for i = 1:length(probs)
%     product = - probs(i) * log2(probs(i));
%     if isnan(product) == 0
%         ntrpy  = ntrpy + product;
%     end
% end
