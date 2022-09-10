# completely algebraic 
processPoly := proc(curve, x, y)
    local m,mm;  
    if indets( evalf(curve),'name') minus {x,y}<>{}  then	#if other symbols are used we got a problem	
		clearCanvas();  
        error "Curve should be a polynomial in both variables %1 and no other variables", {x,y}
    elif  not irreduc(curve) then 
		clearCanvas();  
        error "curve is irreducible"
	elif indets(curve,float)<>{} then #if its a float we got a problem	
		clearCanvas();  
        error "No floating point coefficients allowed"
	fi;
    if member('bug', _Env_GlobalOptions) then
        print("processing polynomial");
    fi;  
    numSheets := degree( curve, y);
    g:=algcurves[genus](curve,x,y);	
	if numSheets = 0 then
		clearCanvas();  
        error "curve is independent of the second variable"
	elif degree( curve, x )=0 then
		clearCanvas();  
        error "curve is independent of the first variable"
	elif g=0 and (member('int',_Env_GlobalOptions) or member('secondary',_Env_GlobalOptions)) then
        clearCanvas();
        error "genus 0 is trivial";
    fi;
    F[1] :=collect( primpart( curve,{x,y}),{x,y}, 'distributed', Normalizer);
    F["y"] := diff(F[1],y); 
    F[2] := evala( normal( resultant(F[1], F["y"], y)  *  lcoeff(F[1],y) ));      
    F[2] := evala( quo( F[2] , gcd( F[2], diff(F[2],x) ), x) );
    F[2] := convert(evala(F[2]/lcoeff(F[2],x)),horner,[x,y]); 
    F[1] := convert(F[1], horner, [x,y]); 
    F["x"] := convert(diff(F[1],x), horner, [x,y]); 
    F["y"] := convert(diff(F[1],y), horner, [x,y]);
    #F["xy"] := diff(F["x"],y); 
    #F["yy"] := diff(F["y"],y);
    #F["xx"] := diff(F["x"],xx);  
    F["Yx"] := -F["x"]/F["y"];   
    #F["Yxx"] := convert(
    #    ( -(F["xx"]+F["xy"]*F["Yx"])*F["y"]+ 
    #F["x"]*(F["xy"]+ F["yy"]*F["Yx"] ))
    #/F["y"]^2, horner, [x,y]);   
     
    #m roots reflected around center. dist can be no less than mm, angle made is half of a sector of angle 2*pi/m
    # so angle dist reps how close we are to center 
    
   
end proc:
getOrderedRoots := proc(x,y,p,x0)
    local approx;
    Digits:=UI_HDigits;
    approx:=toFloat~(convert(polySolve( subs(x=x0,p),y),list));
    if member('bugOrder', _Env_GlobalOptions ) then
        print("bugOrder");
        print(cat("Approximate y at x=", convert(x0, string),
            ". With Digits=", convert(Digits,string)));  
        print(approx);            
    fi;    
    approx:=sortLex(approx);
    if member('bugAC3', _Env_GlobalOptions ) then
        print(cat("After sort: " ));  
        print(approx);                   
    fi;
    return approx;
end proc;

NewtonForRational:=proc(p,q,py,qy,y,y0, DI)
    local P:=p,Q:=q, Py:=py, Qy:=qy,yOld:=y0,yNew,mm:=infinity,num,denum;
    Digits:=minDigits;
    while mm<10^(-DI) do
        if Digits<DI then
            Digits:=2*Digits;
        fi;
        num:= evalf(subs(y=yOld, convert(y*(Py*Q-P*Qy)-P*Q,horner,y)));
        denum:= evalf(subs(y=yOld, convert( Py*Q-P*Qy,horner,y)));    
        yNew:=evalf(num/denum);
        mm:=evalf(ABS(yNew,yOld));
        yOld:=yNew;
    end do;
    return yNew;
end proc;

newtonOnCurve:=proc(p,y,y0)
    local i,II,py:=diff(p,y),yOld:=y0,yNew,mm:=infinity,DI;
    Digits:=UI_HDigits;
    if member('bugNewton',_Env_GlobalOptions) then
        II:=1;
        print("~~~NewtonCheck");
        print("y0=",y0);
    fi;
    while mm>10^(-UI_Digits-EXTRA1) do               
        yNew:=[seq(subs(y=i,y-p/py),i=yOld)]; 
        print(yNew);
        print(yOld);
        mm:=max(abs~(yNew-yOld));
        print(yNew[1]-yOld[1]);
        yOld:=yNew;
        if member('bugNY',_Env_GlobalOptions) then
            II:=1+II;
            print(cat(convert(II,string),": ",convert(mm,string)));
        fi; 
    end do;
    if member('bugNY',_Env_GlobalOptions) then
        print(cat("Took ", convert(II,string), " itterations",
        " to approximate"));
    fi;
    return convert(yNew,list);
end proc;
newtonY:=proc(p,y,y0)
    local II,py:=diff(p,y),yOld:=y0,yNew,mm:=infinity,DI;
    DI:=Digits;
    Digits:=2;
    if member('bugNY',_Env_GlobalOptions) then
        II:=1;
    fi;
    while mm>10^(-DI) do
        if member('bugNY',_Env_GlobalOptions) then
            II:=1+II;
        fi;
        if Digits<DI+2 then
            Digits:=2*Digits;
        fi;        
        yNew:=subs(y=yOld,y-p/py); 
        mm:= evalf(abs(yNew-yOld));
        yOld:=yNew;
    end do;
    return yOld;
end proc;

#given p(z), an itial guess list, and accuracy D get all roots with accuracy D
DK := proc( p, z::name, initalGuess::list,cm,fibre)
    local i,j, roots0,roots1, m:=nops(initalGuess),P, 
	I_m:={seq(i,i=1..m)}, d0,d1,W, w;  
    Digits:=Digits+EXTRA1;
    if lcoeff(p) <>1 then
        P:=evalf(p/lcoeff(p));
    else
        P:=evalf(p);  
    fi;
    roots0 := initalGuess;                
	d0 := mindist( complexLToPtL(roots0) )[1];
    W :=  [seq( evalf( subs( z = roots0[i],
             P) / product( roots0[i]-roots0[ (I_m minus {i})[j] ] , j = 1..m-1 )
        ), i=1..m ) ];
    w := max({seq(ABS(W[i],0),i=1..m)});      
    if w >= cm*d0 then
        if member('bugAC3',_Env_GlobalOptions )  then 
		    print("~~~AC3");
		    print("FAIL ");
	        print("the fibs had min dist ",d0);
            print("w val,", w);
            print("polynomial,", P);   
        fi; 	
        return FAIL;
    fi; 
    roots1 := evalf( [seq(roots0[i]-ABS(W[i],0), i=1..m)] ); #gets you one better guess
	d1:=mindist(complexLToPtL(roots1) )[1];
     if member('bugDK',_Env_GlobalOptions )  then 
		print("~~~AC3");
		print("roots where approximated ");
	    print("they had min dist ",d0);
        print("after new ittereate this dist is,", d1);
        print("ratio", evalf(d1/d0));    
      
        print("want fibre,",fibre);
        print("dont want roots,",roots1);
        print(p);
        print(fsolve(p,z,complex));
        print(P);
        print(fsolve(P,z,complex));
       
    fi;
    rootsOfP:=permuteList(fibre,roots1); #gives back the accurate fibre at that pt 
 
    return rootsOfP,d1/d0;
end;
