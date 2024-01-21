function y = differential(x,n)
for i=1:length(x)-n;
    y(i)=x(i+n)-x(i);

end
end

