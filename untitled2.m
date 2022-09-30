% MATLAB code to plot circle using centre and radius.
  
% create a white image of size 300X600
I=zeros(300, 600)+1;
  
% Radius of circle
R=50;
  
% coordinates of centre point
c=300;
r=150;
  
% accessing every pixel
for i=1:size(I, 1)
    for j=1:size(I, 2)
        dist=round(sqrt((j-c)^2+(i-r)^2));
        if dist==R
            I(i, j)=0;
        end
    end
end
  
% display the image
figure
plot(dist,);