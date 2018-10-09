% Example of using CaptureFigVid
% Cheers, Dr. Alan Jennings, Research assistant professor, 
% Department of Aeronautics and Astronautics, Air Force Institute of Technology

%% Set up 3D plot to record
figure(171);clf;
surf(peaks,'EdgeColor','none','FaceColor','interp','FaceLighting','phong')
daspect([1,1,.3]);axis tight;

%% Set up recording parameters (optional), and record
OptionZ.FrameRate = 20;
OptionZ.Duration  = 15;
OptionZ.Periodic  = true;
% CaptureFigVid([-120,30;-20,10;-110,10;-190,80;-290,10;-380,10], 'WellMadeVid',OptionZ)
% CaptureFigVid([-120,30;-80,30;-40,30;0,30; 40,30;80,10], 'Shelf-DS',OptionZ)
% CaptureFigVid([14,8;54,8;94,8;134,8; 174,8;214,8; 254,8; 294,8; 334,8; 360,8], 'Shelf-Segm',OptionZ)
CaptureFigVid([-110,18;-80,18;-40,18;0,18; 40,18;80,18], 'Branding-Segm',OptionZ)