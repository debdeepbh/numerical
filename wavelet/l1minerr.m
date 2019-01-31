clear all
testvec_6;

% don't forget to run cvx_setup in the installation path before running
% this codez

testvec = [testvec, zeros(1,512)];
%suppL = 10; %sum(abs(testvec) > 0.1)



n = 1024;
fM = dftM(n);

for samplecount = 1:5
%samplecount = 1;
    sampleS(samplecount) = 5+5*samplecount;
    
   %figure (1);
    %    hold on; 
    for count=1:5
    subset = randperm(1024, sampleS(samplecount));
    fMat = fM(subset,:);
        
        cvx_begin
              variable x(n);
             % total variation norm
            minimize( norm(x - circshift(x, [1,0]),1) );
            %minimize( norm(x,2) );
            subject to
                fMat*x == fMat*(testvec.');
        cvx_end
        err(count) = norm(x - testvec',1); %sum((x - testvec') > 0.001);
        xout(:,count) = x;
        
        %plot(x)
        
    end
    incorr(samplecount) = mean(err);
    %xlim([1 n]);
    %hold off
    
%figure (2);
%hold on
%plot(testvec);

 %   plot(mean(xout'))
  %  xlim([1 n]);
    
   % hold off
    %legend('original','average reconstructed')
%    suppL
end
plot(sampleS, incorr)