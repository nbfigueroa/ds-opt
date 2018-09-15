function [Data, Data_sh, att, x0_all, dt] = processDrawnData(data)
M = size(data{1},1)/2;
Data = []; x0_all = []; x0_end = []; Data_sh = [];
for l=1:length(data)    
    data_ = data{l};
    x0_end = [x0_end data_(1:M,end)];
    Data = [Data data_];
    x0_all = [x0_all data_(1:M,1)];
    
    % Shift data to origin for (O2)
    data_(1:M,:) = data_(1:M,:) - repmat(data_(1:M,end), [1 length(data_)]);
    data_(M+1:end,end) = zeros(M,1);

    Data_sh = [Data_sh data_];
end

% Global Attractor of DS
att = mean(x0_end,2);

% dt of Dataset
data_12 = data{1}(:,1:M);
dt = abs((data_12(1,1) - data_12(1,2))/data_12(3,1));
end