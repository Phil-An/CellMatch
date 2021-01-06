function D = euclideanD(A,B)
% USE ONLY ROW VECTORS!
sA = size(A,1);
sB = size(B,1);

if sA == sB
    for i=1:sA
        D(i,:) = sqrt((A(i,1)-B(i,1))^2 + (A(i,2)-B(i,2))^2);
    end
end

if sA == 1 && sA < sB
    for i=1:sB
        D(i,:) = sqrt((A(1,1)-B(i,1))^2 + (A(1,2)-B(i,2))^2);
    end
end

if sB == 1 && sA > sB
    for i=1:size(sA)
        D(i,:) = sqrt((A(i,1)-B(1,1))^2 + (A(i,2)-B(1,2))^2);
    end
end

end