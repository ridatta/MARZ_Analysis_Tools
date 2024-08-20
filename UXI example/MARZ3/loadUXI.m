function out = loadUXI(inDir,camera,frame)

    fname = dir(fullfile(inDir,['*PHC' num2str(camera) '*_09*_F' num2str(frame) '.tif']));
    fname_base = dir(fullfile(inDir,['*PHC' num2str(camera) '*Baseline_F' num2str(frame) '.tif'])); % filenames
    
    img = imread([inDir, fname.name]); % load images
    img_base = imread([inDir, fname_base.name]);
    sub_img = img - img_base; % subtact
    
    if camera == 2
       sub_img = rot90(sub_img,2);
       sub_img = fliplr(sub_img);
    end
    
    % zero
     a = sub_img(:,1:20);
     sub_img = sub_img - mean(a(:));
    
    out = sub_img; 
end