function [image1 image2]= matchimagesize(image1, image2)


% get image dimensions
i1=size(image1,1);
i2=size(image2,1);
i12=size(image1,2);
i22=size(image2,2);

% get max image dimensions
m1=max(i1,i2);
m2=max(i12,i22);

% check for color images
image1rgb=0;
image2rgb=0;

if size(image1,3)== 3
    image1rgb=true;
end

if size(image2,3)== 3
    image2rgb=true;
end



if (i1 < m1)
    d=m1-i1;
    if image1rgb
        image1 = [zeros(ceil(d/2),i12,3); image1; zeros(floor(d/2),i12,3)];
    else
        image1 = [zeros(ceil(d/2),i12); image1; zeros(floor(d/2),i12)];
    end
end

if (i2 < m1)
    d=m1-i2;
    if image2rgb
        image2 = [zeros(ceil(d/2),i22,3); image2 ; zeros(floor(d/2),i22,3)];
    else
        image2 = [zeros(ceil(d/2),i22); image2 ; zeros(floor(d/2),i22)];
    end
    
end

if (i12 < m2)
    d=m2-i12;
    if image1rgb
        image1 = [permute(zeros(ceil(d/2),m1,3),[2 1 3]) image1 permute(zeros(floor(d/2),m1,3),[2 1 3])];
    else
        image1 = [zeros(ceil(d/2),m1)' image1 zeros(floor(d/2),m1)'];
    end
end

if (i22 < m2)
    d=m2-i22;
    if image2rgb
        image2 = [permute(zeros(ceil(d/2),m1,3),[2 1 3]) image2 permute(zeros(floor(d/2),m1,3),[2 1 3])];
    else
        image2 = [zeros(ceil(d/2),m1)' image2 zeros(floor(d/2),m1)'];
    end
end

end
