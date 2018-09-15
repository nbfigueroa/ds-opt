function [Data, Data_sh, att, x0_all, data, dt] = load_dataset_DS(pkg_dir, dataset, sub_sample, nb_trajectories)

dataset_name = [];
switch dataset
    case 1
        dataset_name = '2D_messy-snake.mat';
    case 2
        dataset_name = '2D_Lshape.mat';
    case 3
        dataset_name = '2D_Ashape.mat';
    case 4
        dataset_name = '2D_Sshape.mat';
    case 5
        dataset_name = '2D_multi-behavior.mat';
    case 6
        dataset_name = '3D_viapoint_2.mat';
        traj_ids = [1 2];
    case 7
        dataset_name = '3D_sink.mat';
    case 8
        dataset_name = '3D_bumpy-snake.mat';
    case 9 
        dataset_name = '3D_Cshape_top.mat';
    case 10
        dataset_name = '3D_Cshape_bottom.mat';               
end

if isempty(sub_sample)
   sub_sample = 2; 
end

% For the messy-snake dataset which is already at the origin
if dataset == 1
    Data_ = load(strcat(pkg_dir,'/datasets/',dataset_name));
    data = Data_.data;
    Data = Data_.Data(:,1:sub_sample:end);
    Data_sh  = Data;
    x0_all = Data_.x0_all;
    att = [0 0]';
    data_12 = data{1}(:,1:2);
    dt = abs((data_12(1,1) - data_12(1,2))/data_12(3,1));

% Processing for the 2D Datasets
elseif dataset <= 5
    data_ = load(strcat(pkg_dir,'/datasets/',dataset_name));
    data = data_.data;    
    % Computing the attractor and shifting all the trajectories
    N = length(data);att_ = [];
    M = size(data{1},1)/2;
    for n=1:N   
        att_ = [att_ data{n}(1:M,end)];
    end
    att = mean(att_,2);
    shifts = att_ - att;
    Data = []; Data_sh = []; x0_all = [];
    for l=1:N
        % Gather Data
        data_ = data{l}(:,1:sub_sample:end);
        shifts_ = repmat(shifts(:,l),[1 length(data_)]);
        data_(1:M,:)       = data_(1:M,:) - shifts_;
        data_(M+1:end,end) = zeros(M,1);
        data_(M+1:end,end-1) = (data_(M+1:end,end) + zeros(M,1))/2;
        data_(M+1:end,end-2) = (data_(M+1:end,end-2) + data_(M+1:end,end-1))/2;
        Data = [Data data_];
        
        % All starting position for reproduction accuracy comparison
        x0_all = [x0_all data_(1:2,1)];    
        
        % Shift data to origin for Sina's approach + SEDS
        data_(1:2,:) = data_(1:2,:) - repmat(att, [1 length(data_)]);
        data_(3:4,end) = zeros(2,1);
        
        Data_sh = [Data_sh data_];
        
        % Generate new data structure for SEDS + Diff-DS
        data{l} = data_;
    end
    data_12 = data{1}(:,1:2);
    dt = abs((data_12(1,1) - data_12(1,2))/data_12(3,1));
    
% Processing for the 3D Datasets
else
    data_ = load(strcat(pkg_dir,'/datasets/',dataset_name));
    data = data_.data;
    N = length(data);    
    Data = [];
    for l=1:length(traj_ids)
        % Gather Data
        data_ = data{traj_ids(l)};
        Data = [Data data_(:,1:sub_sample:end)];
    end
    dt = 0.01;
end
end