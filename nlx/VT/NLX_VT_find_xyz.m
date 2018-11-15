function [ x1, y1, z1] = find_xyz( x_a, y_a, x_b,  y_b);
% Function: find_xyz
% input: the location vectors in 2Dfrom the two cameras - x_a y_a from camera 1 and x_b y_b from camera 2
% output: the location of the object in 3D 
% description: the function uses the DLT algorithem to do the transofrmation from 2 vectors of 2D into 1 vector of 3D 
% if the vectors from the cameras are not in the same size it returns an error
% if one of the values in the 2D vectors is 0, it suggests the 3D location will be NaN
%  the 3d x, y and z.
% All DLT-related calculation use a coordinate system with (0,0,0) at top
% left corner of the room --> the final coordinated are transformed
% relative to bottom-left corner of the room (on the floor). This transformation was slightly corrected on 29/3/2011.

% on 28.02.2011 a transformation was added from the measured space to the
% fixed one, using two_tforms.mat and applying on line 53.

load two_tforms; % <- 28.02.2011 updated line

size1 = size(x_a,2);
size2 = size(y_a,2);
size3 = size(x_b,2);
size4 = size(y_b,2);
if not(isequal(size1,size2,size3,size4))
    error ('The input vectors from the cameras are not with the same size');
end

x1 = zeros(1,size1);
y1 = zeros(1,size1);
z1 = zeros(1,size1);

%load mega-calibration data = "fix_array2":
load fix_array2  %  Columns 4:7 of fix_array2 = Replaced for each camera x-->x' and y--> y' corrected for lens distortions 
Cut_pos1=[]; % Empty array means "use all mega-calib points"
F_pos = fix_array(1:246,1:3);
L_pos_c1=fix_array(1:246,4:5); % Camera 1 mega-calib data

% Use control points to calculate the coefficients needed for DLT
% ("learning phase")
 [A_pos_c1,avgres_pos_c1] = DLTFU (F_pos,L_pos_c1,Cut_pos1);
 L_pos_c2=fix_array(1:246,6:7); % Camera 2 mega-calib data
 [A_pos_c2,avgres_pos_c2] = DLTFU (F_pos,L_pos_c2,Cut_pos1);
  
  A_pos=[A_pos_c1,A_pos_c2]; %Coeffs based on Camera 1,2
 
  

 for i=1:size1   % Loop over all frames -- check whetehr there are 0's (one camera did not see) -- convert 0's to NaN's
     % another thing we do is to change the coordinates orientation.
     % coordinates are given relative to the foam coated walls
     % the zero point of the system is in the close left bottom corner of the
     % room, the positive direction is up (Z,col3) , and to the right (X,col1)
     % and towards the inside of the room (Y,col2)

     L=[tforminv(tform1, x_a(i),y_a(i)) tforminv(tform2, x_b(i),y_b(i))]; % <- 28.02.2011 updated line
%    L=[x_a(i), y_a(i), x_b(i),  y_b(i)]; % <- the original line of the previous one

    if not(x_a(i)&x_b(i)&y_a(i)&y_b(i))
        x1(i) =NaN;
        y1(i)= NaN;
        z1(i) = NaN;
    else
        ans1 = RECONFU(A_pos,L); %Use cooff. calculated above to get 3D coordinates
               
        x1(i) = ans1(1); % the x axis in the room is with good cooridnates - no need to fix,
        y1(i)= ans1(2)+4600; % this is the y axis in the room  [29/3/2011 -- shift corrected from 4590 to 4600 mm - Nachum+Michael]
        z1(i) = ans1(3)+2660; % this is the z axis in the room  [29/3/2011 -- shift corrected from 2585 to 2660 mm - Nachum+Michael]
        
    end
    
 end
  
 