function [ ] = PSD_Circles( Data, Centers, Labels, varargin)
%function [ ] = Dist_Line( Dist, Labels )
%   Detailed explanation goes here

% [x,y] = ginput(23);
% c=[x y];
% c=[Xcenter Ycenter];
if length(varargin) == 1
    Visibility = varargin{1};
else
    Visibility = 'on';
end

r=10;
% figure('Visible', 'off','units','normalized','outerposition',[0 0 1 1])
figure('units','normalized','outerposition',[0 0 1 1], 'Visible', Visibility)
% set(gca,'ydir','reverse');
circles(120,[0 0],'color','y','linewidth',2)
circles(70,[0 0],'color','y','linewidth',2)
circles(40,[0 0],'color','y','linewidth',2)
circles(20,[0 0],'color','y','linewidth',2)

circles(r,Centers,'color','b','linewidth',0.5)
axis equal
hold on

for i=1:length(Labels)
    text(Centers(i,1),Centers(i,2),Labels{i},'HorizontalAlignment','center','FontSize',10);
end

m = mean(Data);

for i = 1:length(Data)
    if Data(i) > m
        for j = 1:10
            circles(r+j-5,Centers(i,:),'color','k','linewidth',0.5);
        end
    end
end

% count=0;
% for i=1:length(Labels)
%     for j=i+1:length(Labels)
%         if Dist(i,j) < 0.2 
%             count=count+1;
%             line([Center(i,1) Center(j,1)],[Center(i,2) Center(j,2)],'color',[5*(.2-Dist(i,j)) 5*(.2-Dist(i,j)) 5*(.2-Dist(i,j))],'LineWidth',((1-((Dist(i,j))*5))));
%         end
%     end
% end
% count = count / ((length(Labels)*(length(Labels)-1))/2)

end

