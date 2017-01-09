function [cname cext] =  denoi(filen, band)    % works for color and grayscale
Im = imread(filen);
[path, namepart, ext] = fileparts(filen);

if band == 0
    chan = Im;
else
    chan = Im(:,:,band);
end

disp('Adding noises.');
Noi(:,:,1) = imnoise(chan,'gaussian');
Noi(:,:,2) = imnoise(chan,'poisson');
Noi(:,:,3) = imnoise(chan,'salt & pepper');


for i=1:3

    noisyFileName = strcat(namepart,int2str(band),'-noise-',int2str(i),ext);
    imwrite(Noi(:,:,i),noisyFileName);
    fprintf('File with noise type %d saved to: %s\n', i, noisyFileName);
        fprintf('Applying denoising method');
        % gaussian, susan, ynf, median, NL
        fprintf(' Gaussian filter,');
        DeN(:,:,1) = gaussF2(5,2,Noi(:,:,i));
        fprintf(' SUSAN filter,');
        DeN(:,:,2) = susanF2(5,100,Noi(:,:,i));
        fprintf(' YNF,');
        DeN(:,:,3) = ynf(2,250,Noi(:,:,i));
        fprintf(' Median filter');
        DeN(:,:,4) = medF(2,Noi(:,:,i));
        
        
        szz = size(DeN);
        fprintf('Denoised files saved to: \n');
        for j = 1: szz(3)
            denoiseFile = strcat(namepart,int2str(band),'-denoise-',int2str(i),'-',int2str(j),ext);
            imwrite(DeN(:,:,j),denoiseFile);
            fprintf('%s\n',denoiseFile);
        end
        fprintf('.\n');
end
cname = namepart;
cext = ext;


% denoise




% imshow(ig)