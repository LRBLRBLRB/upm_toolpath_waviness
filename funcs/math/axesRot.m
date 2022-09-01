function R = axesRot(varargin)
% to solve the rotation matrix based on the known coordinate axes
%
% usage:
% R = axesRot(x1,y1,x2,y2,'xy');
% R = axesRot(x1,y1,x2,y2,'xy');

switch varargin{end}
    case 'xy'
        x1 = varargin{1}./norm(varargin{1});
        y1 = varargin{2}./norm(varargin{2});
        z1 = cross(x1,y1);
        z1 = z1./norm(z1);
        x2 = varargin{3}./norm(varargin{3});
        y2 = varargin{4}./norm(varargin{4});
        z2 = cross(x2,y2);
        z2 = z2./norm(z2);
        R = [x2'*x1, y2'*x1, z2'*x1;
            x2'*y1, y2'*y1, z2'*y1;
            x2'*x1, y2'*z1, z2'*z1];
    case 'yz'
        y1 = varargin{1}./norm(varargin{1});
        z1 = varargin{2}./norm(varargin{2});
        x1 = cross(y1,z1);
        x1 = x1./norm(x1);
        y2 = varargin{3}./norm(varargin{3});
        z2 = varargin{4}./norm(varargin{4});
        x2 = cross(y2,z2);
        x2 = x2./norm(x2);
        R = [x2'*x1, y2'*x1, z2'*x1;
            x2'*y1, y2'*y1, z2'*y1;
            x2'*x1, y2'*z1, z2'*z1];
    case 'zx'
        z1 = varargin{1}./norm(varargin{1});
        x1 = varargin{2}./norm(varargin{2});
        y1 = cross(z1,x1);
        y1 = y1./norm(y1);
        z2 = varargin{3}./norm(varargin{3});
        x2 = varargin{4}./norm(varargin{4});
        y2 = cross(z2,x2);
        y2 = y2./norm(y2);
        R = [x2'*x1, y2'*x1, z2'*x1;
            x2'*y1, y2'*y1, z2'*y1;
            x2'*x1, y2'*z1, z2'*z1];
    case 'xyz'
        x1 = varargin{1}./norm(varargin{1});
        y1 = varargin{2}./norm(varargin{2});
        z1 = varargin{3}./norm(varargin{3});
        x2 = varargin{4}./norm(varargin{4});
        y2 = varargin{5}./norm(varargin{5});
        z2 = varargin{6}./norm(varargin{6});
        R = [x2'*x1, y2'*x1, z2'*x1;
            x2'*y1, y2'*y1, z2'*y1;
            x2'*x1, y2'*z1, z2'*z1];
    otherwise
        R1 = vecRot(varargin{1},varargin{3});
        varargin{1} = R1*varargin{1};
        varargin{2} = R1*varargin{2};
        R2 = vecRot(varargin{2},varargin{4});
        R = R2*R1;
end

end