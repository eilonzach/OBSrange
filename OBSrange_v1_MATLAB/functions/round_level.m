function [ xx ] = round_level( x, y )
%[xx] = round_level(x,level)
% rounds the number/vector x to the nearest y (the level)
xx = round(x./y).*y;
end

