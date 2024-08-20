function addScale(scale,img,xtickvals,ytickvals)
% scale = px/mm
[m,n] = size(img);
% convert to px
xtickvals_px = xtickvals / scale + n /2; 
ytickvals_px = ytickvals / scale + m /2; 
xticks(round(xtickvals_px));
yticks(round(ytickvals_px));
xticklabels(xtickvals);
yticklabels(ytickvals);
end