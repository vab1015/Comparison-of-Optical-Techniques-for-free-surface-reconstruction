%% Computing free surface elevations using Optic flow and TIR-D
%  See Jain et al. (2021) for details of TIR-D setup and equations used herein. 
%  See https://www.mathworks.com/help/vision/ref/opticalflowfarneback.html
%  for info about MATLAB's version of Optic flow- Farneback method

%  invgrad2_rect and designgrad1D by Wildeman (2018) is used for least-square integration. Retrieved from: 
%https://www.mathworks.com/matlabcentral/fileexchange/67258-fast-checkerboard-demodulation-for-digital-schlieren-imaging

%  Prepared by Vivek Bheeroo. 
%  Last modified on 01/27/2023
clear
clc
%% Loading video and computing displacement field 
tic
vidReader = VideoReader('TIRD_image_stack_W1.avi');  %Video is an image stack with every odd numbered frame being the reference checkerboard. 
opticFlow= opticalFlowFarneback; %Setting optic flow to Farneback method

%Tuning the parameters of Farneback. See detailed description at : 
%https://www.mathworks.com/help/vision/ref/opticalflowfarneback.html
opticFlow.NumPyramidLevels = 3;  %Number of pyramid levels. 
                                 %Default =3, range= [1,inf)
opticFlow.PyramidScale = 0.5;  %Rate of downsampling at each pyramid level. 
                               %Default = 0.5, range= (0,1) 
opticFlow.NumIterations = 3;   %Number of search iterations. 
                               %Default= 3, range- [1,inf)
opticFlow.NeighborhoodSize= 10;  %Size of the pixel neighborhood. 
                                 %Default= 5, range- [1,inf)
scale_filt = 3;                  %Multiplier for filter size (next line)

opticFlow.FilterSize = floor(15*scale_filt);  %Filter size. Default = 15. 
%range= [2,inf). Multiply by scale_filt to move up by multiples of 15

%Running function 'estimateFlow' to compute Farneback algorithm
frameGray = im2gray(readFrame(vidReader)); %Initializing gray scale matrix 
flowp = estimateFlow(opticFlow,frameGray); %Initializing optical flow output 
while hasFrame(vidReader)
    frameRGB= readFrame(vidReader);
    frameGray = rgb2gray(frameRGB);  
    flowp(end+1) = estimateFlow(opticFlow,frameGray);
end 

Num_frames =vidReader.NumFrames; 
N = Num_frames-1; %Output from estimateFlow will be N-1 cells
U_cell_data = cell(N,1); 
V_cell_data = cell(N,1);
%Loop is run in a way so that only even elements of flowp is considered.
%This is because the odd numbered frames are the reference image. 
for jj=1:2:N
    flowss = flowp(1,jj+1);
    vxx = flowss.Vx;
    vyy = flowss.Vy;
    U_cell_data{jj,1} = flipud(vxx);      %To match MATLAB convention 
    V_cell_data{jj,1} = flipud(vyy.*-1);  %To match MATLAB convention 
end
final_U_cell_data = U_cell_data; 
final_V_cell_data = V_cell_data; 
final_U_cell_data = final_U_cell_data(~cellfun('isempty',final_U_cell_data)); %Gets rid of empty cells
final_V_cell_data = final_V_cell_data(~cellfun('isempty',final_V_cell_data)); %Gets rid of empty cells 
N = length(final_U_cell_data); 

clearvars U_cell_data; 
clearvars V_cell_data; 

L_x = 0.104; %Streamwise domain size [m]
L_y = 0.065; %Transverse domain size [m]
%Constructing x and y arrays and meshgrid 
X = linspace(0,L_x,size(final_U_cell_data{1,1},1)); 
Y = linspace(0,L_y,size(final_V_cell_data{1,1},2)); 
[x,y] = meshgrid(X,Y);

%% Computing surface gradients and integrating to obtain final surface elevations
theta= 62; %Incident angle, according to Jain et al.(2021) 

%Modiying displacement field. See Eqn 15 from Jain et al.(2021)
u_f = final_U_cell_data; 
v_f = cellfun(@(x) x.*(cosd(theta))^2,final_V_cell_data,'un',0); 
clearvars final_U_cell_data
clearvars final_V_cell_data

H= 0.175; %Water depth [m]
T= 0.5*sind(2*theta); %Eqn 15

fact = exp(y.*T./H); %Multipling factor used in of RHS of eqn 21

dx = x(1,2) - x(1,1);  %resolution in x [m]
dy = y(2,1) - y(1,1);  %resolution in y [m]

t= linspace(0,(N/50),N); %Creating a time vector 
eta = cell(N,1);  %Initializing surface elevation cell array 
for jj=1:N
    u_calib = u_f{jj,1}.*(L_x/size(u_f{1,1},1)); %Converting x displacements from pixel to real-world coordinates
    v_calib = v_f{jj,1}.*(L_x/size(u_f{1,1},1)); %Converting y displacements from pixel to real-world coordinates
    fx = u_calib .*-(fact'/(2*H));  %RHS of Eqn 21 (see Jain et al. (2021)) 
    fy = v_calib .*-(fact'/(2*H));  %RHS of Eqn 21 (see Jain et al. (2021)) 
    fhat = invgrad2_rect(-fy,-fx);  %Numerical integration. Negatives sign to match convention used in Wildeman (2018)
    eta{jj,1} = (fhat .*dx)./fact'; 
    sprintf('Tracking progress of numerical integration: %d/%d complete',jj,N) 
end 
computation_time = toc

%% Creating animation of free surface elevations - Tune caxis, zlim and campos based on Exp
myVideo = VideoWriter('TIRD Free-surface');
myVideo.FrameRate= 5;
open(myVideo)
figure('units','normalized','outerposition',[0 0 1 1])
for jj=1:100
    colormap('viridis') 
    hh=surf(x,y,(flipud(eta{jj,1}')));
    set(hh,'EdgeColor','none')
    set(gca,'FontSize',25)
    title(sprintf('%d/%d frames, %.2f s',jj,N,t(1,jj)))
    aa=colorbar;
    ylabel(aa,'\eta [m]','FontSize',16)
    caxis([-2*10^-5 2*10^-5])
    zlim([-3*10^-5 3*10^-5])
    campos([-0.3802   -0.5582    0.0013])
    pbaspect([x(1,end)/y(end,1) 1 1])

    drawnow
    pause(0.0001)

    frame = getframe(gcf);
    writeVideo(myVideo, frame);
end
close(myVideo)

%% Saving displacement field (in pixels) and surface elevations 
TIRD_output.U_disp = u_f; 
TIRD_output.V_disp = v_f; 
TIRD_output.eta = eta; 

save TIRD_output.mat TIRD_output -v7.3 
