function [] = GIFmaker( gcf, framenum, Directory, Name, varargin)
%function [] = GIFmaker( gcf, framenum, Directory)
%function [] = GIFmaker( gcf, framenum, Directory, DelayTime)
%   Detailed explanation goes here
    if length(varargin) == 1
        Delay = varargin{1};
    else
        Delay = 0.6;
    end
    filename = [Directory Name '.gif'];
    frame = getframe(gcf);
     im = frame2im(frame);
     [A,map] = rgb2ind(im,256); 

     if framenum == 1
		imwrite(A,map,filename,'gif','LoopCount',1,'DelayTime',Delay);
	 else
		imwrite(A,map,filename,'gif','WriteMode','append','DisposalMethod','leaveInPlace','DelayTime',Delay);
     end

end

