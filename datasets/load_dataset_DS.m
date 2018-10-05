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
        dataset_name = '3D_viapoint_3.mat';
    case 7
        dataset_name = '3D_sink.mat';
    case 8 
        dataset_name = '3D_Cshape_bottom.mat';
    case 9
        dataset_name = '3D_Cshape_top.mat';                       
    case 10
        dataset_name = '3D-pick-box.mat';                       
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
    N = length(data);
    for l=1:N
        % Gather Data
        data{l} = data{l}(:,1:sub_sample:end);
    end
    [Data, Data_sh, att, x0_all, dt, data] = processDataStructure(data);
    
% Processing for the 3D Datasets
else
    data_ = load(strcat(pkg_dir,'/datasets/',dataset_name));
    data_ = data_.data;
    N = length(data_);    
    data = []; dt = 0.01;
    traj = randsample(N, nb_trajectories)'
    for l=1:nb_trajectories
        % Gather Data
        data{l} = data_{traj(l)}(:,1:sub_sample:end);
    end
    [Data, Data_sh, att, x0_all, ~, data] = processDataStructure(data);
end
end