function handleCol(filen)   % feed image file

Im = imread(filen);
si = size(Im);
szise = size(si);

if szise(2) == 2    % it has only one layer
    [namepart ext] = denoi(filen,0);
else
    for band = 1:si(3)
        fprintf('Layer=%d\n',band);
        [namepart ext] = denoi(filen,band);
    end
    
    % Construct color noisy
    disp('Constructing color images...');
    for i=1:3
        for band=1:si(3)
            Noi(:,:,band) = imread(strcat(namepart,int2str(band),'-noise-',int2str(i),ext));
        end
        noiName = strcat(namepart,'-col-',int2str(i),ext);
        imwrite(Noi,noiName);
        fprintf('Color image with noise type %d saved to: %s', i, noiName);
        for j=1:4
            for band=1:si(3)
                Noi(:,:,band) = imread(strcat(namepart,int2str(band),'-denoise-',int2str(i),'-',int2str(j),ext));
            end
            fprintf(' %d,',j);
            colOut = strcat(namepart,'-col-',int2str(i),'-denoise-',int2str(j),ext);
            imwrite(Noi,colOut);
            fprintf('Denoising type %d applied to noise type %d saved to: %s',j,i,colOut);
        end
        fprintf('\n');
    end
    

end




