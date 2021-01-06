function cropped3Darray = imcrop3D (array, rec)
%% crop 3D Matrix imcrop3D(stack, rec)
recx=rec(3);
recy=rec(4);


if and(recx >= size(array(:,:,1),2), recy >= size(array(:,:,1),1))
    tI = zeros(recy,recx );
    for i=1:size(array,3)
        [array_matched(:,:,i) ~]= matchimagesize(array(:,:,i), tI);
    end
      
elseif recx >= size(array(:,:,1),2)
    d = recx -size(array(:,:,1),2);
    d= d + recx + floor(recx*0.1);
    tI = zeros(recy, d);
    
    for i=1:size(array,3)
        [array_matched(:,:,i) ~]= matchimagesize(array(:,:,i), tI);
    end
    
elseif recy >= size(array(:,:,1),1)
    d = recy -size(array(:,:,1),1);
    d= d + recy + floor(recy*0.1);
    tI = zeros(d, recx);  
     
    for i=1:size(array,3)
        [array_matched(:,:,i) ~]= matchimagesize(array(:,:,i), tI);
    end
else
    array_matched = array;
end


for i=1:size(array,3)
    cropped3Darray(:,:,i)=imcrop(array_matched(:,:,i), rec);
end

end