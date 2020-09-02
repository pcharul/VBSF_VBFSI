function [err] = mre_error(true,predicted,mask)

el=find(~true);
    mask(el)=1;
    mask = ~mask;
    true1=true.*mask;
    predicted1=predicted.*mask;
    pos_tst=find(mask);
    %%

    err=norm(true1(pos_tst)-predicted1(pos_tst),'fro')/norm(true1(pos_tst),'fro');
end
