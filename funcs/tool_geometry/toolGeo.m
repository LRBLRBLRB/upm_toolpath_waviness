% 刀具几何(toolGeo)类
classdef toolGeo
    % properties 存放类的属性
    properties
        center;
        radius;
        profileMeas;
    end
    % methods 存放类的成员函数
    methods
        function item = toolGeo(Center,Radius,ProfileMeas)
            if nargin == 2
                item.profileMeas = zeros(0,0);
            else
                item.profileMeas = ProfileMeas;
            end
            item.center = Center;
            item.radius = Radius;
        end
    end
end