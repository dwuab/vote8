function Pt=three_state_evolve( W12, p, eta, N, P0 )

Pt=zeros([N+1 N+1]);
Pt2=zeros([N+1 N+1]);
dt=1;
Pt(P0(1)*N+1,P0(2)*N+1)=1;
T=3000;

F(T/10) = struct('cdata',[],'colormap',[]);

for t=1:T
    x=1; y=1; z=N+1;
    Pt2(x,y)=Pt(x,y)-dt*(W12(x,y)+W12(x,z)+W12(y,x)+W12(y,z) ...
        +W12(z,x)+W12(z,y))*Pt(x,y);
    
    Pt2(x,y)=Pt2(x,y)+dt*W12(x+1,z-1)*Pt(x+1,y)+dt*W12(y+1,z-1)*Pt(x,y+1);
    for i=2:N
        x=i; z=N-(y-1)-(x-1)+1;
        Pt2(x,y)=Pt(x,y)-dt*(W12(x,y)+W12(x,z)+W12(y,x)+W12(y,z) ...
            +W12(z,x)+W12(z,y))*Pt(x,y);
        Pt2(x,y)=Pt2(x,y)+dt*W12(x+1,z-1)*Pt(x+1,y)+dt*W12(y+1,x-1)*Pt(x-1,y+1)...
            +dt*W12(y+1,z-1)*Pt(x,y+1)+dt*W12(z+1,x-1)*Pt(x-1,y);
    end
    x=N+1; z=N-(y-1)-(x-1)+1;
    Pt2(x,y)=Pt(x,y)-dt*(W12(x,y)+W12(x,z)+W12(y,x)+W12(y,z) ...
        +W12(z,x)+W12(z,y))*Pt(x,y);
    Pt2(x,y)=Pt2(x,y)+dt*W12(y+1,x-1)*Pt(x-1,y+1)+dt*W12(z+1,x-1)*Pt(x-1,y);
    
    for j=2:N
        y=j; x=1; z=N-(y-1)-(x-1)+1;
        Pt2(x,y)=Pt(x,y)-dt*(W12(x,y)+W12(x,z)+W12(y,x)+W12(y,z) ...
            +W12(z,x)+W12(z,y))*Pt(x,y);
        Pt2(x,y)=Pt2(x,y)+dt*W12(x+1,y-1)*Pt(x+1,y-1)+dt*W12(x+1,z-1)*Pt(x+1,y)+...
            dt*W12(y+1,z-1)*Pt(x,y+1)+dt*W12(z+1,y-1)*Pt(x,y-1);
        for i=2:N-j+1
            x=i; z=N-(y-1)-(x-1)+1;
            Pt2(x,y)=Pt(x,y)-dt*(W12(x,y)+W12(x,z)+W12(y,x)+W12(y,z) ...
                +W12(z,x)+W12(z,y))*Pt(x,y);
            Pt2(x,y)=Pt2(x,y)+dt*W12(x+1,y-1)*Pt(x+1,y-1)+dt*W12(x+1,z-1)*Pt(x+1,y)...
                +dt*W12(y+1,x-1)*Pt(x-1,y+1)+dt*W12(y+1,z-1)*Pt(x,y+1)+...
                dt*W12(z+1,x-1)*Pt(x-1,y)+dt*W12(z+1,y-1)*Pt(x,y-1);
        end
        x=N-j+2; z=N-(y-1)-(x-1)+1;
        Pt2(x,y)=Pt(x,y)-dt*(W12(x,y)+W12(x,z)+W12(y,x)+W12(y,z) ...
                +W12(z,x)+W12(z,y))*Pt(x,y);
        Pt2(x,y)=Pt2(x,y)+dt*W12(x+1,y-1)*Pt(x+1,y-1)+dt*W12(y+1,x-1)*Pt(x-1,y+1)...
            +dt*W12(z+1,x-1)*Pt(x-1,y)+dt*W12(z+1,y-1)*Pt(x,y-1);
    end
    
    x=1; y=N+1; z=1;
    Pt2(x,y)=Pt(x,y)-dt*(W12(x,y)+W12(x,z)+W12(y,x)+W12(y,z) ...
        +W12(z,x)+W12(z,y))*Pt(x,y);
    Pt2(x,y)=Pt2(x,y)+dt*W12(x+1,y-1)*Pt(x+1,y-1)+dt*W12(z+1,y-1)*Pt(x,y-1);
    
    Pt=Pt2;
    
    if rem(t,10)==0
    surf(Pt,'EdgeColor','none');
    view([0 0 1]);
    
    F(t/10)=getframe;
    end
end

%movie(F);

end

