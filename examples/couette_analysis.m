load ../output/data/couette_010_3_1.mat

t_all = t;
rot_all = rot;
xc_all = xc;
tau_all = tau;

load ../output/data/couette_010_3_1_tmp.mat
t_all = [t_all, t_all(end) + t(2:end)];
rot_all = [rot_all, rot];
tau_all = [tau_all;tau(2:end,:)];
xc_all = cat(3, xc_all, xc(:,:,2:end));

%load ../output/data/couette_010_3_1_part3.mat
%t_all = [t_all, t_all(end) + t(2:end)];
%rot_all = [rot_all, rot];
%tau_all = [tau_all;tau(2:end,:)];
%xc_all = cat(3, xc_all, xc(:,:,2:end));

%load ../output/data/couette_010_3_1_part4.mat
%t_all = [t_all, t_all(end) + t(2:end)];
%rot_all = [rot_all, rot];
%tau_all = [tau_all;tau(2:end,:)];
%xc_all = cat(3, xc_all, xc(:,:,2:end));

% compute angles relative to radial direction
angles = zeros(prams.np, length(t_all));
r = zeros(prams.np,length(t_all));

for i = 1:length(t_all)
   r_angles = atan2(xc_all(2,:,i),xc_all(1,:,i));
   
   angles(:,i) = r_angles - wrapToPi(tau_all(i,:));
   
   r(:,i) = sqrt(xc_all(2,:,i).^2 + xc_all(1,:,i).^2) - 5;
end


plot(t_all(1:end-1)/pi,rot_all/837.7580)
xlabel('percentage of outer cylinder revolution');
title('Torque on inner cylinder');


%ylim([0,max(rot)]);

% h = figure();
% addpath('../src/matlab2tikz/src');       
%              
% for it = 1:length(t_all(t_all < pi))
%     
%     if it == 1
%         mkdir('../output/gifs/angle_hist_010');
%     end
%     
%     y = radtodeg(mod(wrapToPi(angles(:,it)),pi));
%     
%     histogram(y,20,'Normalization','probability');
%     
%     title(sprintf('percentage of outer revolution = %6.3f', t_all(it)/pi));
%     xlabel('Angle');
%     ylabel('Percentage of particles');
%     ylim([0,0.4]);
%     xlim([0,180]);
%     
%     drawnow;
%     pause(0.1);
%     frame = getframe(h);
%     im = frame2im(frame);
%     [imind, cm] = rgb2ind(im,256);
%     
%     cleanfigure();
%     matlab2tikz(['../output/gifs/angle_hist_010/couette_010_angles', sprintf('%04d', it), '.tikz'],...
%         'height', '10cm', 'width', '12cm', 'standalone', true, 'floatFormat', '%.3f');
%     
% end
% hold on
% y = abs(mod(radtodeg(wrapToPi(angles')),180) - 90);
% z = zeros(size(t_all));
% 
% for i = 1:size(y,2)
%     col = y(:,i)';  % This is the color, vary with x in this case.
%     surface([t_all;t_all]/pi,[y(:,i),y(:,i)]',[z;z],[col;col],...
%             'facecol','no',...
%             'edgecol','interp',...
%             'linew',2);
% end
% 
%     
% % hold on
% % colMap = parula(256);
% % [np, n_t_steps] = size(angles);
% % for i = 1:n_t_steps
% %     for j = 1:np
% %         
% %         dotColor = colMap(ceil(abs(mod(radtodeg(wrapToPi(angles(j,i)+1e-9)),180)-90)*256/90),:);
% %         plot(t_all(i)/pi, abs(mod(radtodeg(wrapToPi(angles(j,i))),180)), '.', 'Color', dotColor);
% %     end
% % end
% %plot(t_all/pi, mod(radtodeg(wrapToPi(angles')),180), 'b.');
% xlabel('percentage of outer cylinder revolution');
% title('Particle angles');
% 
% figure()
% hold on
% y = r';
% z = zeros(size(t_all));
% 
% for i = 1:size(y,2)
%     col = y(:,i)';  % This is the color, vary with x in this case.
%     surface([t_all;t_all]/pi,[y(:,i),y(:,i) ]',[z;z],[col;col],...
%             'facecol','no',...
%             'edgecol','interp',...
%             'linew',2);
% end
% % hold on
% % % colMap = parula(256);
% % [np, n_t_steps] = size(angles);
% % for i = 1:n_t_steps
% %     for j = 1:np
% %         
% %         dotColor = colMap(ceil((r(j,i))*256/5),:);
% %         plot(t_all(i)/pi, r(j,i), '.', 'Color', dotColor);
% %     end
% % end
% %plot(t_all/pi, r', 'b');
% xlabel('percentage of outer cylinder revolution');
% title('Particle distance from inner cylinder');
