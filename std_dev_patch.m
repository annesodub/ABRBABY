function [] = std_dev_patch(grpA, grpB, std_err, color)


curve1 = grpA.
% y = time vector; d = data vector; C = color

% y = [1,2,3,4,5];
% d = [1,3.4,5,7.5,0.5];
% e= [31, 12, 0.1, 6.9, 15];
% C = 'r';

%y = rand(1,10); % your mean vector;
x = 1:numel(data);
%stderror= std(data) / sqrt( length( data ))
%std_dev = 1;
%std_dev =b;
curve1 = d + std_dev;
curve2 = d - std_dev;
x2 = [x, fliplr(x)];
inBetween = [curve1, fliplr(curve2)];
%fill(x2, inBetween, 'Color', color, 'FaceAlpha', 0.1, 'EdgeColor', 'none');
fill(x2, inBetween, C, 'FaceAlpha', 0.1, 'EdgeColor', 'none');
%patch(curve1, fliplr(curve2), 'g');
hold on;
%plot(x, d, 'Color', color, 'LineWidth', 2);
%plot(x, d, C, 'LineWidth', 2);


%Script from function HaloPatchMAD
%mad = median(abs(d-repmat(median(d'),size(d,2),1)')');
%mad = std(d)/sqrt(length(d));

% Set mean, low and high halo borders
% avg = mean(d,2) ;
% Lhi = mean(d,2) + mad' ; 
% Llow = mean(d,2) - mad' ;

% Draw halo patch 
% hPatch = patch([vTime fliplr(vTime)], [Lhi'  fliplr(Llow')], [0.6 0.6 0.6],...
%         'FaceAlpha', 0.1, ...
%         'FaceColor', C/255, ...
%         'EdgeColor', 'none', ...
%         'Parent',    hAxes);
% 
% hPatch = patch([vTime fliplr(vTime)], [Lhi'  fliplr(Llow')], [0.6 0.6 0.6],...
%         'FaceAlpha', 0.1, ...
%         'FaceColor', C, ...
%         'EdgeColor', 'none');
% 
%     % Skip the name of the previous plot from the legend
% hPatch.Annotation.LegendInformation.IconDisplayStyle = 'off';



%figure; plot(y,d,C);hold on; plot(y, e, 'b');hold on ; std_dev_patch(y, d, C);hold on ; std_dev_patch(y, e, 'b')
end 
