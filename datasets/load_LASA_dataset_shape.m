function [demos, name] = load_LASA_dataset_shape(number_model)
% This is a matlab function illustrating 30 human handwriting motions
% recorded from Tablet-PC. These datas can be found in the folder
% 'DataSet'. 
%
% Please acknowledge the authors in any academic publications
% that have made use of this library by citing the following paper:
%
%  S. M. Khansari-Zadeh and A. Billard, "Learning Stable Non-Linear Dynamical 
%  Systems with Gaussian Mixture Models", IEEE Transaction on Robotics, 2011.
%
% To get latest upadate of the software please visit
%                          http://lasa.epfl.ch/khansari
%
% Please send your feedbacks or questions to:
%                           mohammad.khansari_at_epfl.ch


%%
names = {'Angle','BendedLine','CShape','DoubleBendedLine','GShape',...
         'heee','JShape','JShape_2','Khamesh','Leaf_1',...
         'Leaf_2','Line','LShape','NShape','PShape',...
         'RShape','Saeghe','Sharpc','Sine','Snake',...
         'Spoon','Sshape','Trapezoid','Worm','WShape','Zshape',...
         'Multi_Models_1', 'Multi_Models_2', 'Multi_Models_3','Multi_Models_4'};
     
name = names{number_model};

if number_model<0 || number_model>30
    disp('Wrong model number!')
    disp('Please try again and type a number between 1-30.')
elseif number_model == 0
    return
else
    %% preprocessing
    load(['DataSet/' names{number_model}],'demos','dt') %loading the model
end

end