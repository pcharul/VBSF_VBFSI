function [Y] = add_outlier_mtx(Y,s)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

load('datasets/outlier_delhi.mat');
[m,n]=size(Y);
dat=find(Y);
dat=dat(randperm(length(dat)));
Obs = dat(1:round(s*length(dat)));
ind = Obs;
sz = [m,n];
[row,col] = ind2sub(sz,ind);
Y(ind)=Y(ind)+outlier(row);


end

