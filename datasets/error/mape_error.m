function [err] = mape_error(true,predicted,mask)
el=find(~true);
    mask(el)=1;
    mask = ~mask;
    true1=true.*mask;
    predicted1=predicted.*mask;
    pos_tst=find(mask);
    
err=sum(abs(true1(pos_tst)-predicted1(pos_tst)))./length(pos_tst);
  
end

