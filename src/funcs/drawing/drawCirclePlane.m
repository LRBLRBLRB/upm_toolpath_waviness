function varargout = drawCirclePlane(varargin)
%drawPlane: To draw a plane


center = varargin{1};
radius = 1.2*varargin{2};
normVec = varargin{3};
normVec = normVec/norm(normVec);
plotOpts = varargin{4};

xLim = [center(1) - radius, center(1) + radius];
yLim = [center(2) - radius, center(2) + radius];

X = [xLim(1),xLim(2),xLim(2),xLim(1)];
Y = [yLim(1),yLim(1),yLim(2),yLim(2)];
Z = zeros(1,4);
Z(1) = center(3) - ...
    (normVec(1)*(X(1) - center(1)) + normVec(2)*(Y(1) - center(2)))/normVec(3);
Z(2) = center(3) - ...
    (normVec(1)*(X(2) - center(1)) + normVec(2)*(Y(2) - center(2)))/normVec(3);
Z(3) = center(3) - ...
    (normVec(1)*(X(3) - center(1)) + normVec(2)*(Y(3) - center(2)))/normVec(3);
Z(4) = center(3) - ...
    (normVec(1)*(X(4) - center(1)) + normVec(2)*(Y(4) - center(2)))/normVec(3);

varargout{1} = fill3(X,Y,Z,plotOpts.FaceColor, ...
    'FaceAlpha',plotOpts.FaceAlpha,'EdgeColor',plotOpts.EdgeColor);

end