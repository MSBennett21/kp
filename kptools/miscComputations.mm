 
 # This procedure gives the permutation that transforms 
# lijst1 to lijst2. Because of numerical roundoff, equalities
# between elements are never used.

getPermutation:=proc(lijst1::list,lijst2::list)
        local outlijst, m,mm,ii, jj,r,S,d;
        Digits:=UI_HDigits;
        if member('bugPerm',_Env_GlobalOptions) then
            print("permuting", lijst1);
            print("into", lijst2);
        fi; 
        outlijst:=NULL;
        S:={seq(ii,ii=1..nops(lijst1))};
        for ii in lijst2 do
            m:=infinity;
            for jj in S minus {outlijst} do
                mm:=ABS(ii,lijst1[jj]);
                if evalf(mm-m)<0 then
                    d:=m;
                    m:=mm;
                    r:=jj;
                elif evalf(mm-d) <0 then
                    d:=mm;
                fi;
            od;
            outlijst:=outlijst,r;
            #if fnormal(evalf(m*10-d)) >0 then
            #    print(m*10-d);
            #    error "Computation inaccurate, use larger value for Digits"
            #fi;
        od;
        #print(outlijst);
        [outlijst]
    end;
permuteList:=proc(lijst1,lijst2)
    local outlijst, m,mm,ii,j, jj,r,S,d,l1;
    outlijst:=[];
    l1:={op(lijst1)};
    for ii from 1 to nops(lijst2) do
        m:=infinity;
        for jj in l1 do
            mm:=ABS(lijst2[ii],jj);
            if evalf(mm-m)<0 then               
                m:=mm;
                j:=jj;                
            fi;
        od; 
        l1:=l1 minus{j};
        outlijst:=[op(outlijst),j];
    od;
    outlijst;
end;
permuteFibs:=proc(perm,lijst::list)
    local i;
    [seq(lijst[perm[i]],i=1..nops(lijst))];
end;
wIsToRightOfz:= proc( z0, w0  ) 
    local Z0:=convert(z0,rational), W0:=convert(w0,rational); 
    W0:=normal(W0*(Re(Z0)-I*Im(Z0))/(Re(Z0)^2+Im(Z0)^2));
    Z0:=normal(Z0*(Re(Z0)-I*Im(Z0))/(Re(Z0)^2+Im(Z0)^2));
    
    if member('bugCross',_Env_GlobalOptions ) then 
        print("~~~~CROSS");
        print(cat("the primary direction ", convert(z0,string)));
        print(cat("the secondary direction ", convert(w0,string)));
        print(cat("primary scale",convert(Z0,string))); 
        print(cat("secondary scale",convert(W0,string))); 
        print(plots[display]({
            plot({[[0,0],toFloat(complexToPt(z0))],[[0,0],toFloat(complexToPt(w0))]},color="red"),
            plot({[[0,0],toFloat(complexToPt(Z0))],[[0,0],toFloat(complexToPt(W0))]},color="blue")},
            axes=none,title="old is red"));                  
        print(normal(Im( W0 )));
    fi;
    return evalb(normal(Im( W0 ))< 0);
end proc;    
nPtBoundary := proc( ll::list,n:=6,  scale:=2)
    local center,k,r,N:=nops(ll); 
    Digits:=UI_HDigits;
    try
	    center :=  scrub(add( ll[k], k=1..N )/N);
    catch:
        error "given list is not of the right form: must be a list of complex numbers";
    end try;
    
    if not ( type(n, integer) and scale > 0 )then
        error "wrong form of input"
    fi;
	r :=  scale*maximalDistWRT(ll,center); 
     
    boundaryPoints:=table();
    boundaryPoints["radius"]:= r;
    boundaryPoints["center"]:= center;
    boundaryPoints["boundary"] :=  sortArg(
        scrub~( [ seq(
        center + r * exp( (2/n*Pi*k+Pi/n)*I),
        k=1..n )]), -r); 
end proc:
maximalDistWRT := proc( ll::list,  c::complexcons)
   return max(  map( a ->  abs(  a- c ), ll))  ;
end proc:
scrubToList := proc()
    local i;        
    if nargs> 0 then
        try
    	    return 
            [ seq( 
                ifelse( type(args[i],indexable),
                    scrub~(args[i]),scrub(args[i])),
                i=1..nargs) ];
            catch:                
                print(args);
                error "scrubing the given data did not work";
            end try;
    else       
        error "no input given";  
    fi;    
end:
isNormalForm :=proc(e) ::boolean;
    try 
        return evalb(e[1]=edgeData[{op(e)}]["normDir"][1]);
    catch:
        print(e);
        error "wrong use of normal form"
    end try;
end proc; 
ABS:=proc(z1,z2)
    return toFloat(sqrt((Re(z1)-Re(z2))^2+(Im(z1)-Im(z2))^2));
end proc;
clearCanvas := proc()
    #initalized data
    infDist:=1;
    minStep:=1/2^BaseLvl;
    BLvl:=BaseLvl;
    basePoint:=infinity; 
    basePreimage:=infinity; 
    numOfProblemPts:=0;
    GN:=0;GLvl:=0;
    minDigits:=0; plotter:={}; F:=table(); 
    edgeData:=table(); vertexData:=table();
    pathData:=table(); E:={};V:={};
    UI_Digits:=Digits;UI_HDigits:=Digits;UI_LDigits:=Digits;
    unassign('homBasis', 'problemPoints',
        'singPoints', 'boundaryPoints', 'branchPoints', 'diffBasis',  'g', 'numSheets',
        'dxx','dxxx', 'intData','homBasis', 'dxx','dxxx','intData');   
   
      	
    gc();
end:
toFloat := proc(val)
    Digits:=Digits+EXTRA1;
	return simplify( fnormal( evalf(val),Digits-EXTRA1), zero);
end:
scrub := proc(val)
    local mm;
    return toFloat(floor( toFloat(10^(Digits)*val))/10^(Digits));
end:
complexToPt := proc( z )       
    return ifelse( z <> infinity, [Re(z),Im(z)],NULL);
end:
complexLToPtL := proc( ll:: list)
    local i;        
    return [seq(complexToPt(ll[i]),i=1..nops(ll))];
end:
ptLToComplexL := proc( ll)
    return map(a->complex(a[1],a[2]),ll);
end:
constructEdge := proc( ll, b )
	return map(a->{b,a}, ll);
end:
edgeToPrint := proc( ll )
	return map(a->complexLToPtL([a[1],a[2]]), ll);
end:
sortLex:=proc(ll)
    local func;  
    Digits:=UI_HDigits;
	func := (s::complexcons,t::complexcons) -> evalb( 
             Re(s) < Re(t) or 
            ( Re(s) = Re(t) and ( Im(s) <= Im(t))) 
            
            );
	sort(ll, func )
end;
sortArg:=proc(ll,b::complexcons)    
    #need the b that is the focal point of this ordering
    Digits:=UI_HDigits;
	func :=(s::complexcons,t::complexcons) -> evalb( 
        s=b or 
        (t <> b and argument( s - b) < argument( t - b )) or 
        ( argument( s-b ) = argument( t-b ) and  ABS( s,b ) <= ABS( t , b ) ) ) ; 
	sort(ll,func )	
end;

lineIntersection:=proc(z1,z2,w1,w2)    
   local  x1:=Re(z1), y1:=Im(z1), x2:=Re(z2), y2:=Im(z2), 
   x3:=Re(w1),  y3:=Im(w1), x4:=Re(w2),  y4:=Im(w2),xx,yy,D;
    
    D:=toFloat((x1-x2)*(y3-y4)-(y1-y2)*(x3-x4));     
    xx:=((x1*y2-x2*y1)*(x3-x4)-(x1-x2)*(x3*y4-y3*x4));
    yy:=((x1*y2-x2*y1)*(y3-y4)-(y1-y2)*(x3*y4-y3*x4));
    if member('bugV3',_Env_GlobalOptions ) then 
        print("denominator of intersection pt",D);
        print(toFloat((xx+I*yy)/D));
    fi;  
    return toFloat((xx+I*yy)/D);
end;
updatePlotter := proc()
    local listOfPointSets,listOfEdgeSets,ptSet,eSet,e;
    if nargs=0 then
        error "no input into plotter";
    elif nargs=1 then
        listOfPointSets:=args[1]; 
        listOfEdgeSets:=[];
    else

        listOfPointSets:=args[1];
        listOfEdgeSets:=args[2];
    fi;  
    if nops(args[1][1])>0 then  
        for ptSet in listOfPointSets do
            plotter := plotter union { plots[pointplot](
            complexLToPtL(convert(ptSet[1],list)),op(ptSet[2]))};
        od; 
    fi;
    for eSet in listOfEdgeSets do #{[e1,e2,..],options}
        plotter:= plotter union { plot( op~(eSet[1]),op(eSet[2]))};
    od; 
end;
 getPsuedoSite:=proc(siteIndex)
        return ifelse( siteIndex <= numOfProblemPts-1, problemPoints[siteIndex],boundaryPoints["boundary"][siteIndex-numOfProblemPts+1]);
    end;
getSite:=proc(siteIndex)
	    return ifelse( siteIndex <= numOfProblemPts-1, problemPoints[siteIndex],infinity);
    end;
getSiteIndex:=proc(site)
    local i,II:=0,m:=infinity;
    Digits:=UI_HDigits;
        if site=infinity then
            return numOfProblemPts;
        else 
            for i from 1 to numOfProblemPts-1 do
                if ABS(problemPoints[i],site)<m then
                    II:=i;
                    m:=ABS(problemPoints[i],site);
                fi;
            od;
	        return II;
        fi;
    end;
