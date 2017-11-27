function [skypos] = imagepos_to_skypos(platescale,pa,x_offset,y_offset)
%IMAGEPOS_TO_SKYPOS - convert image position in x and y pixel offset from center to RA and Dec offset from center on sky 
%
% Syntax:  [skypos] = imagepos_to_skypos(platescale,pa,x_offset,y_offset)
%
% Inputs:
%   platescale - Instrument plate scale in mas/pix
%   pa - Position Angle of detector array in degrees measured counter-clockwise from y axis to true direction of north
%   x_offset - x offset from center of image in pixels
%   y_offset - y offset from center of image in pixels
%
% Outputs:
%   skypos - [dRA dDec] Position on sky in milliarcsecond (mas) offset from center
%
% Example: 
%   sky_position = imagepos_to_skypos(plate_scale,PA,dx,dy)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files or other data files required: none
%
% See also: N/A

% Author: Ricky Nilsson, Ph.D., Astronomy
% Department of Astronomy, California Institute of Technology, CA, USA
% email address: ricky@caltech.edu  
% Website: http://www.caltech.edu/~rnilsson
% November 2017; Last revision: 27-Nov-2017

%------------- BEGIN CODE --------------

pos_offset = [x_offset y_offset 0];
t = (pa/180)*pi; % position angle in radians

% Rotation matrix for plane rotation by angle t around z axis
Rz = [cos(t) -sin(t) 0; sin(t) cos(t) 0; 0 0 1];

% Scaling matrix
S = [1 0 0; 0 1 0; 0 0 0] * platescale;

S_pos_offset = pos_offset * S; % Scaled position offset
RzS_pos_offset = S_pos_offset * Rz; % Scaled and rotated position offset

skypos = RzS_pos_offset(1:2);


%------------- END OF CODE --------------