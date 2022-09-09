
#FACT if 0 in one direction functionis zero function
checkWCError := proc( RM::Matrix, U::list, V::list, W::list, d::complexcons, c := 0, DI:=2 )
    local mm, i, bool:=false;
   # if eq(Re(RM), zero( nops(U) ))  then
   #     for i from 1 to nops(U) do
    #        if fnormal(evalf((U[i]-conjugate(U[i]))/2), Digits-6) <> 0.0 then 
    #            bool := true;
    #            break;
    #        fi;
    #    od;        
    #fi;
    
    mm:=-1;
    if member('line',{args}) then
        mm:= abs(HBF(0, 0,0, RM, U, V, W, d, c, DI ));
        #print(mm);
        for i from 1 to 20 do
            if abs(HBF( (i-1)/19, 0, 0, RM, U, V, W, d, c, DI )) >= mm then 
                mm := abs(HBF((i-1)/19, 0,0, RM, U, V, W, d, c, DI ));
                #print(mm);
            fi;            
        od;
    else
        mm:=max([seq(abs( HBF( evalf(rand()/(10^12 - 1)), evalf(rand()/(10^12 - 1)), evalf(rand()/(10^12 - 1)), RM, U, V, W, d, c, DI )), i=1..20)] ); 
    fi;
    mm;    
end proc;
HBF := proc(x,y,t, RM::Matrix, U::list, V::list, W::list, d::complexcons, c:=0, DI:=2 )
    local argument := convert( evalf(<op(U)> * x + <op(V)> * y + <op(W)> * t), list),
    Theta, ThetaX, ThetaT, ThetaY, 
    ThetaXX, ThetaXT, ThetaYY, 
    ThetaXXX, ThetaXXXX, tol := 0.1^DI; 
    #print(argument);
    #print(RM, [], tol, output = list);
    Theta := RiemannTheta( argument, RM, [], tol, output = list )[2];
    ThetaX := RiemannTheta( argument, RM, evalf([  U  ]), tol, output = list )[2];
    ThetaY := RiemannTheta( argument, RM, evalf([ V ]), tol, output = list )[2];
    ThetaT := RiemannTheta( argument, RM, evalf([ W ]), tol, output = list )[2];
    ThetaXX := RiemannTheta( argument, RM, evalf([ U, U ]), tol, output = list )[2]; 
    ThetaXT := RiemannTheta( argument, RM, evalf([ U, W ]), tol, output = list )[2];
    ThetaYY := RiemannTheta( argument, RM, evalf([ V, V ]), tol, output = list )[2];
    ThetaXXX := RiemannTheta( argument, RM, evalf([ U, U, U ]), tol, output = list )[2];
    ThetaXXXX := RiemannTheta( argument, RM, evalf([ U, U, U, U ]), tol, output = list )[2];

   
    evalf[DI](ThetaXXXX * Theta - 
    4 * ThetaXXX * ThetaX + 
    3 * ThetaXX^2 + 
    4 * ThetaX * ThetaT - 
    4 * ThetaXT * Theta + 
    3 * ThetaYY * Theta - 
    3 * ThetaY^2 + 
    8 * d * Theta^2 +
    6*c*ThetaXX*Theta -
    6*c*ThetaX^2);      
     
end proc;


