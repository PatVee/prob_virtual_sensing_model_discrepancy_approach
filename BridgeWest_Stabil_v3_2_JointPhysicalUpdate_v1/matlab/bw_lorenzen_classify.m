function T=bw_lorenzen_classify(fe)
% A screening classification, not an experimental mode identification.
nm=numel(fe.freq);index=(1:nm)';freq=fe.freq(:);family=strings(nm,1);
lobes=zeros(nm,1);symmetric=lobes;vertical=lobes;topRatio=lobes;
mask=fe.x>fe.c.bearing_x_m(1)+1e-6 & fe.x<fe.c.bearing_x_m(2)-1e-6;
for j=1:nm
    A=fe.P(6*(fe.top(1,mask)-1)+3,j);B=fe.P(6*(fe.top(2,mask)-1)+3,j);
    energy=A'*A+B'*B;symmetric(j)=sum((A+B).^2)/(2*energy+1e-30);
    U=reshape(fe.P(:,j),6,[])';vertical(j)=sum(U(:,3).^2)/sum(sum(U(:,1:3).^2));
    topRatio(j)=max(abs([A;B]))/(max(abs(U(:,3)))+1e-30);
    if symmetric(j)>=.5,v=A+B;else,v=A-B;end
    v=v(abs(v)>.08*max(abs(v)));lobes(j)=1+sum(v(1:end-1).*v(2:end)<0);
    if vertical(j)>.5 && topRatio(j)>.25
        if symmetric(j)>=.5,family(j)="bending";else,family(j)="torsion";end
    else,family(j)="local/lateral";end
end
T=table(index,freq,family,lobes,symmetric,vertical,topRatio,'VariableNames', ...
    {'fe_index','frequency_Hz','family','lobes','symmetric_fraction','vertical_fraction','top_amplitude_fraction'});
end
