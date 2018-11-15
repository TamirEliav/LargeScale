function side = side_of_line(B1,B2,points)
    d = (points(:,1)-B1(1)).*(B2(2)-B1(2)) - (points(:,2)-B1(2)).*(B2(1)-B1(1));
    side = sign(d);
end