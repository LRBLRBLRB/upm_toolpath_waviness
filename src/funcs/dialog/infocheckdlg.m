function isContinue = infocheckdlg(workspaceDir,toolFileName,unit,aimRes,arcLength,maxAngPtDist,rMax,c)
%INFOCHECKDLG 

fig = uifigure('Name','Surface Generation Infomation Check', ...
                'WindowStyle','alwaysontop','Visible','off');
fig.CloseRequestFcn = @(app,event)fig_close_req(app);
fig.Position(3:4) = [400,300];
figGridLayout = uigridlayout(fig,[2,2]);
figGridLayout.RowHeight = {'1x','fit'};
figGridLayout.ColumnWidth = {'1x','1x'};
color = fig.Color;

txt = uitextarea('Parent',figGridLayout,...
   'Value',{'Surface was generated successfully!','', ...
    'The parameters are listed below:', ...
    sprintf('0. Workspace: %s',getlastfoldername(workspaceDir)), ...
    sprintf('1. Tool file: %s',toolFileName), ...
    sprintf('2. Aimed residual error: %f %s',aimRes,unit), ...
    sprintf('3. C increment arc %f %s',arcLength,unit), ...
    sprintf(['   max-angle %f',char(176),')'],maxAngPtDist*180/pi), ...
    sprintf('4. Surface radius: %f %s',rMax,unit), ...
    sprintf('5. Surface curvature: %f%s^{-1}',c,unit),'', ...
    'Ready to continue?'},'WordWrap','on','Editable','off', ...
    'BackgroundColor',color);
txt.Layout.Row = 1;
txt.Layout.Column = [1,2];

figEnterButton = uibutton(figGridLayout,'push','Text','OK & Continue','Visible','on');
figEnterButton.Layout.Row = 2;
figEnterButton.Layout.Column = 1;
figEnterButton.ButtonPushedFcn = @(app,event)fig_enter_button_pushed;

figCancelButton = uibutton(figGridLayout,'push','Text','Cancel','Visible','on');
figCancelButton.Layout.Row = 2;
figCancelButton.Layout.Column = 2;
figCancelButton.ButtonPushedFcn = @(app,event)fig_refit_button_pushed;

fig.Visible = 'on';

uiwait(fig);

    function fig_close_req(app)
        isContinue = 0;
        delete(app);
    end

    function fig_enter_button_pushed
        isContinue = 1;
        fig_close_req(fig);
    end

    function fig_refit_button_pushed
        isContinue = 0;
        fig_close_req(fig);
    end
end