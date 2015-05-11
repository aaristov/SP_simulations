function crop = cropimage(image, x,y,cropsize)
crop=imcrop(image,[x-cropsize/2,y-cropsize/2,cropsize,cropsize]);
end