function month=spotmonth(day)
    for id= 1 : length(day)
        monthDay=[31,28,31,30,31,30,31,31,30,31,30,31];
        monthFirstDay(1)=1;
        for i = 2 : 12
           monthFirstDay(i)=sum(monthDay(1:i-1))+1; 
        end
        temp=day(id)./monthFirstDay;
        month(id)=max(find(temp>=1));
    end
end