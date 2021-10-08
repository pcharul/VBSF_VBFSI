function [Y] = add_outlier_mtx(Y,s)


[m,n]=size(Y);
dat=find(Y);
dat=dat(randperm(length(dat)));
Obs = dat(1:round(s*length(dat)));
ind = Obs;
sz = [m,n];
[row,col] = ind2sub(sz,ind);
m=randi([-100 100],length(row),1);
Y(ind)=Y(ind)+m;


end

