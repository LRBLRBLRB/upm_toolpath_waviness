function projPt = ptOnLine(Pt,varargin)
% to figure out the projection of the given point on the given straight
% line
%
% usage:
%
% projPt = projectionOnLine(Pt,param)
%   the line is described with parameters of the general equation
%   Pt          (3,:)   the point that will be projected
%   param       (3,1)   the parameters of the equation, i.e., A, B & C,
%                       respectively
%   projPt      (3,:)   the projection of the point Pt
%
% projPt = projectionOnLine(Pt,linePt,lineDir)
%   the line is described with a point on it and its direction vector
%   Pt          (3,:)   the point that will be projected
%   linePt      (3,1)   the point on the line
%   lineDir  (3,1)   the direction vector of the line
%   projPt      (3,:)   the projection of the point Pt

switch nargin
    case 1
        param = varargin{1};
        error('it hasn''t been realized');
    case 2
        linePt = varargin{1};
        lineDir = varargin{2};
        ptNum = size(Pt,2);
        if ptNum == 1
            projPt = linePt ...
                + dot(Pt - linePt,lineDir)./norm(lineDir).*lineDir;
        else
            projPt = ndgrid(linePt,1:ptNum) ...
                + dot(Pt - ndgrid(linePt,1:ptNum),ndgrid(lineDir,1:ptNum),1) ...
                ./norm(lineDir).*lineDir;
        end
    otherwise
        error('invalid input format. Only 2 or 3 inputs are allowed.');
end

end