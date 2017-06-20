function res = im2row(im, winSize)
% function to convert 2D subtensors of a tensor to row vectors
% window-by-window
%
% res = im2row(im, winSize)
%
%   Inputs:
%           im      - 4D data array with size [sx, sy, N, D] where
%                     N is the number of phase-cycles, D is the coil
%				      number and the images are sx-by-sy
%           winSize - the size of the window that slides over 2D data
%
%   output:
%           res 	- 4D tensor whose first dimension is the generated
%					  rows, the second dimension is number of windows,
%					  and the other two dimensions remain the same
%
% (based on im2row function of SPIRiT library)

	[sx,sy,N,D] = size(im);

	res = zeros((sx-winSize(1)+1)*(sy-winSize(2)+1),prod(winSize),N,D);
	count=0;
	for y=1:winSize(2)
		for x=1:winSize(1)
			count = count+1;
			res(:,count,:,:) = reshape(im(x:sx-winSize(1)+x,y:sy-winSize(2)+y,:,:),...
				(sx-winSize(1)+1)*(sy-winSize(2)+1),1,N,D);
		end
	end
	
end