function [XL, XR] = LagStructure()
global a c s w0
sr = s;
Xr = c*exp(-(sr/a).^2)+ w0;
XR = [sr', Xr'];

sl = -s;
Xl = c*exp(-(sl/a).^2)+ w0;
XL = [sl', Xl'];


% % %%%%%%%%%==============================================================
% % %%%%%%*************************************************
% global Lx Ly lx0
% YL = zeros(length(s),2);
% YL(:,1) = -s;
% YR = zeros(length(s),2);
% YL(:,1) = s;
% 
% for j = 1:length(s)
%     ss = s(j);
%     if ss <= lx0
% 
%         y = (c+ w0) - (c/lx0)*ss;
% 
%         YL(j,2) = y;
%         YR(j,2) = y;
%     else
%         YL(j,2) = w0;
%         YR(j,2) = w0;
%     end
% end
% 
% figure
% xx = [-Lx, Lx, Lx, -Lx];
% yy = [0, 0, Ly, Ly];
% fill(xx, yy, 'c');
% hold on
% plot(-s, XL(:,2), 'r', s, XR(:,2), 'r','LineWidth',1.5);
% 
% plot(-s, YL(:,2), 'b', s, YR(:,2), 'b','LineWidth',1.5);
% 
% xlim([-Lx, Lx]); ylim([0, Ly]);
% p= gcf;
% exportgraphics(p,'Initial_Geom.png','Resolution',300);

end
