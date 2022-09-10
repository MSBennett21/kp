#fix center problem
#clean up code
# Do very quality test on Hurwitz
#verts--low percision needed
#problem points exact percision needed
#y values - high persision
#geometric high persicion and then scrub to low
Monodromy := proc(curve, x::name, y::name) 
    local mono;
 
    UI_ODigits:= ifelse(
        not assigned(_Env_GlobalOptions) or not member('int',_Env_GlobalOptions),  
        Digits, UI_ODigits);
    UI_Digits := max( 10, Digits);         
    if UI_Digits <> Digits then
        WARNING("Increasing digits to 10");            
    fi;
    _Env_GlobalOptions := ifelse( assigned(_Env_GlobalOptions),
        [ args[ 3..nargs ], op(_Env_GlobalOptions)],[ args[ 3..nargs ]]); 

    processPoly(curve,x,y);
    if member('bugPoly', _Env_GlobalOptions ) then
        print(cat("The function has been processed and is ", convert(F[1], string),
            "." ));   
         print(cat("The polynomial that is procssed is ", convert(F[2], string),
            "." ));            
    fi; 
    mono:= mono1(x,y); 
    
   # if member(group,[args[4..-1]]) then  #if group is mentioned, return the permgroup
	#	C:= procname(args[1..3],op(subs(group=NULL,[args[4..-1]])));
	#	'permgroup'(nops(C[2]),{seq(i[2], i=C[3] )})
	#	[C[1],C[2],C[3]];
   
    if not member('secondary',_Env_GlobalOptions) then 
        clearCanvas();  
    fi; 
    mono;

end:
mono1:=proc(x::name, y::name)   
	if numSheets=1 then
		linearMonodromy(x,y)
	else 
        if degree(F[1], x)< 3 then #come check this out
		    userinfo(1,'algcurves',
		    `Consider switching the roles of the variables. The curve has degree <=2 in the first variable.`);
	    fi;
		generalMonodromy(x,y)
	fi;
end:
# This procedure gives the monodromy information of the Riemann surface
# of a single-valued function.
linearMonodromy:=proc(x::name, y::name)   
    Digits:=UI_Digits;
	basePoint := 0;
	while subs(x = basePoint, lcoeff( F[1], y))=0 do
		basePoint := basePoint - 1;
	od;   
	basePreimage := [  scrub~(fsolve(subs( x=basePoint, curve ), y )) ];
	[basePoint,basePreimage,[]]
end:
generalMonodromy:=proc(x::name, y::name) 
    local Dual,monoData;
	secondInitialization(x,y);
    Digits := UI_Digits;  
    #order of problem points and boundary pts are fixed!
    #should be no more than O(n^2) N  #use the digits that user specified
    if numOfProblemPts>2 then
        Dual := voronoiDual( complexLToPtL( [op(  problemPoints) , op(boundaryPoints["boundary"]) ] ) );      
        produceVoronoi( Dual );
        Digits := UI_HDigits;
        problemPoints := sortArg(problemPoints,basePoint);
        problemPoints := sortArg(problemPoints,basePoint);
        basePreimage :=  scrub~(getOrderedRoots(  x,y,F[1],basePoint));  
        problemPoints := Array(1..numOfProblemPts-1, problemPoints );    
        ArrayTools[Append]( problemPoints, infinity );
        if  member('visual', _Env_GlobalOptions ) or member('showpaths', _Env_GlobalOptions ) or member('bugs', _Env_GlobalOptions ) then
            updatePlotter( [{[basePoint],{symbol=cross, color ="Blue"}}]);                  
            #print(plots[display](op(plotter)));           
        fi;
        if  member('bug', _Env_GlobalOptions ) then
            print(" Vertices, edges, have been produced. Basepoint is determined:",basePoint);
            print("order is now determined on the sheets and on the sites according to lex order");                
            print(basePreimage);
            print("producing hurwitz system");           
        fi;
        produceHurwitz();
        
    else
        problemPoints:=[op(problemPoints),infinity];
        V:=[problemPoints[1]-1,problemPoints[1]-I,problemPoints[1]+1,problemPoints[1]+I];
        basePoint:=V[1];
        for v in V do
            vertexData[ v ] := table();
            vertexData[ v ]["neighbors"] := {op(V)} minus{v};
            vertexData[ v ]["sites"] := problemPoints;
            vertexData[v]["onBoundary"] := true;
        od; 
        
        pathData[problemPoints[1]]:=table();
        pathData[problemPoints[1]]["vCell"]:=V;
        pathData[problemPoints[1]]["path"] :=  [[V[1],V[2]],[V[2],V[3]],[V[3],V[4]],[V[4],V[1]]];
        pathData[problemPoints[2]]:=table();
        pathData[problemPoints[2]]["vCell"]:=V;
        pathData[problemPoints[2]]["path"] :=  [[V[1],V[4]],[V[4],V[3]],[V[3],V[2]],[V[2],V[1]]]; 
        E:={{V[1],V[2]},{V[2],V[3]},{V[3],V[4]},{V[4],V[1]}};
        basePreimage :=  scrub~(getOrderedRoots(  x,y,F[1],basePoint));  
        for i from 1 to 4 do
            edgeData[{op(pathData[problemPoints[1]]["path"][i])}]:=table();
            edgeData[{op(pathData[problemPoints[1]]["path"][i])}]["normDir"]:=pathData[problemPoints[1]]["path"][i];
            edgeData[{op(pathData[problemPoints[1]]["path"][i])}]["gam"] :=normal(convert((pathData[problemPoints[1]]["path"][i][2]-pathData[problemPoints[1]]["path"][i][1])/2,rational));
            edgeData[{op(pathData[problemPoints[1]]["path"][i])}]["path"] :=pathData[problemPoints[1]]["path"][i][1]*(1-s)/2+pathData[problemPoints[1]]["path"][i][2]*(1+s)/2;
        od;
        if  member('visual', _Env_GlobalOptions ) or member('showpaths', _Env_GlobalOptions ) or member('bugs', _Env_GlobalOptions ) then
            updatePlotter( [{[basePoint],{symbol=cross, color ="Green"}}]);                  
            updatePlotter( [{V,{symbol=solidcircle, color = "black"}}]);          
              updatePlotter( 
                    [{}],[{complexLToPtL~(pathData[problemPoints[1]]["path"]), {color="red"}}]);
        fi;
    fi;
    
    
	if member('showpaths', _Env_GlobalOptions ) or member('bug', _Env_GlobalOptions) or  member('visual', _Env_GlobalOptions) then             
        print(plots[display](plotter,axes=NORMAL,scaling=CONSTRAINED, 
		title=typeset("Paths chosen for the analytic continuation of y(x)"), 	
        labels=[Re, Im]));           
    fi;  
    #have egdes, verts, sites
    #no longer need the order on the problem points
    monoData:=generateMonodromy(x,y);    	
	
	
	#	print( plots[display]( plot(  eListToPrint(E) , color ="red" ), 
	#	plots[pointplot]( {op(complexLToPtL(convert(branchPoints,list)))},symbol=solidcircle, color = "black"),
    #    plots[pointplot]( {op(complexLToPtL(convert(singPoints,list)))},symbol=solidcircle, color = "red"),
    #    plots[pointplot]( {op(complexLToPtL(boundaryPoints["boundary"]))},symbol=circle,color = "red" ),
    #    plots[pointplot]( complexToPt(basePoint),symbol=asterisk, color = "red"), );
	#fi;    
	return [basePoint, basePreimage, convert(monoData, list)]; #sites_ordered,P,ACDATA
end proc;
#Input is curve output is f p,q,h
#f =curve
#p = curve that has roots of f that are problem points but not singular
#q = curve that has roots of f that  singular pts
#h = list of sheet cordinates ordered-it induces the order
secondInitialization := proc(x::name, y::name)
    processPoints(x,y);    
    if member('int', _Env_GlobalOptions) then
        intProcessing(x,y);
    fi;   
    if member('bug', _Env_GlobalOptions ) then       
        print(cat("There are ", convert(numOfProblemPts-1, string),
            " finite problem points. We also computed the boundary with center ",
            convert(boundaryPoints["center"], string), " and radius ", convert(boundaryPoints["radius"], string)));
        print( "Working Digits are currently", UI_LDigits,UI_Digits,UI_HDigits);
        print("min dist,",infDist);
        print(cat("everything has been computed to ", convert(UI_HDigits, string)," digits"));
        if member('visual', _Env_GlobalOptions ) then  
            print(plots[display](plotter,axes=none,title="Artificial boundary points and problem points"));
        fi;            
    fi;      	  
end proc;
processPoints:= proc(x,y)
    if member('bug', _Env_GlobalOptions) then
       print("processing points");
    fi;  
#WRITE A CODE THAT TRANSFORMS THE DATA SO THAT IT TIS NICE AND BOUNDED IN A CIRCLE AROUND ZERO
#just use conformal mapings
    Digits := UI_Digits+2*EXTRA1;
    problemPoints := toFloat(convert(polySolve( F[2],x),list));
    numOfProblemPts := nops(problemPoints)+1; 
    #Need to add more generality here
    if numOfProblemPts>2 then
        infDist := mindist(complexLToPtL(problemPoints))[1]/2;
        minDigits :=  max(ceil(-log10(min(infDist,1))),1)+EXTRA2;
        if minDigits  > UI_Digits then 
            clearCanvas();  
            error "Increase digits.";
        elif minDigits> 2*EXTRA1 then
            Digits := UI_Digits+minDigits;
            problemPoints := toFloat(convert(polySolve( F[2],x),list));
        fi;    
        UI_LDigits := minDigits;
        UI_HDigits := UI_Digits + 2*UI_LDigits;  #22 or more
        Digits := UI_HDigits;
        problemPoints := scrub~(problemPoints);       
        nPtBoundary(problemPoints);

        # print(boundaryPoints["radius"]);
        problemPoints :=  sortArg( problemPoints, -boundaryPoints["radius"]);
    
        if  member('showpaths', _Env_GlobalOptions ) or  member('visual', _Env_GlobalOptions ) or member('bug', _Env_GlobalOptions )  then
            updatePlotter( [{problemPoints,{symbol=solidcircle, color = "red"}},
                    {boundaryPoints["boundary"],{symbol=circle,color = "red"}}]);
        fi;
    else
        if  member('showpaths', _Env_GlobalOptions ) or  member('visual', _Env_GlobalOptions ) or member('bug', _Env_GlobalOptions )  then
            updatePlotter( [{problemPoints,{symbol=solidcircle, color = "red"}}]);
           
        fi;
    fi;
end proc;

produceVoronoi := proc(triangleSet)
    local triIndex, sites,realSites,m,mm, v,w, 
    numOfTri:=nops(triangleSet), Ti,
    i,vSites,wSites,j;
   
    if member('bug', _Env_GlobalOptions ) then
        print(cat("There are ", convert(numOfTri,string), " triangles in the Delaunay triangulation"));
        if member('bugV', _Env_GlobalOptions ) then
            print(triangleSet); 
        fi;    
    fi; 
    basePoint:=0;   
    for triIndex from 1 to numOfTri do
		Ti := triangleSet[triIndex]; #get a triangle 
        Digits:=UI_HDigits;
        v := toFloat(intersectionOfBLines( op(getPsuedoSite~(Ti)))); 
       
        if member('bugV', _Env_GlobalOptions ) then
            print("for sites", Ti);
            print("Crt vert hi percision",v);
        fi;
        Digits:=UI_LDigits;
        v:=scrub(v);
        if Re(v)<= Re(basePoint) then
            basePoint:=v;
        fi;
        if member('bugV', _Env_GlobalOptions ) then
            print("low percision",v);
        fi;
        if assigned(vertexData[ v ]) then 
            vertexData[ v ]["sites"][1] := 
            vertexData[ v ]["sites"][1] union getSite~({op(Ti)});
            vertexData[ v ]["sites"][2] := 
            vertexData[ v ]["sites"][2] union getPsuedoSite~({op(Ti)});
            vertexData[v]["onBoundary"] := member(infinity,
            vertexData[ v ]["sites"][1]  );
            if member('bugV', _Env_GlobalOptions )  then  
                print(cat("vertex", convert(v,string),
                " has already been registered. Its current sites are ",
                convert(vertexData[ v ]["sites"][1],string)));
                print(cat("is the vertex on the boundary?", convert(vertexData[v]["onBoundary"],string)));
            fi; 
        else
            vertexData[ v ] := table();
            vertexData[ v ]["neighbors"] := {};
            vertexData[ v ]["sites"] := [getSite~({op(Ti)}), getPsuedoSite~({op(Ti)})];
            vertexData[v]["onBoundary"] := member(infinity,
            vertexData[ v ]["sites"][1]  );
            if member('bugV', _Env_GlobalOptions )  then  
                print(cat("vertex", convert(v,string),
                " has NOT already been registered. Its  sites are ",
                convert(vertexData[ v ]["sites"][1],string)));
                print(cat("is the vertex on the boundary?", convert(vertexData[v]["onBoundary"],string)));
            fi; 
        fi;
        V := V union { v };  
	od;
    infDist:=min(mindist(complexLToPtL(convert(V,list)))[1]/2,infDist);
    if infDist<1/10^(Digits) then  
        print(infDist);
        print(V);
        clearCanvas();  
        error "Increase digits: the verticies of the graph are not being resolved";
    fi;  
    if  member('visual', _Env_GlobalOptions ) 
        or member('showpaths', _Env_GlobalOptions ) then
        updatePlotter( [{convert(V,list),{symbol=solidcircle, color = "black"}}]);
        if (member('bugEdge', _Env_GlobalOptions ) and member('visual', _Env_GlobalOptions ) ) then  
            print(plots[display](plotter,axes=none,title="Verticies have been resolved"));
        fi;            
    fi;    
    mm:=V;
    while nops(mm)> 1 do 
        v:=mm[1];
        vSites:=vertexData[ v ]["sites"]; #should be a tuple of two sets
        if member('bugEdge',_Env_GlobalOptions ) then 
            print(cat("there are ", convert( nops(mm),string)," left"));
            print("considering ", v);
            print(cat("it has the sites ", convert( vSites,string)));
        fi;
        for w in mm[2..-1] do   
            wSites:=vertexData[ w ]["sites"];
            sites:=sortLex(wSites[1] intersect vSites[1]);
            psuedoSites:=wSites[2] intersect vSites[2];
            if nops(psuedoSites)>=2 and nops(sites)>=2 then
                if nops(E) <> nops(E union {{w,v}}) then
                    vertexData[ v ]["neighbors"]:=vertexData[ v ]["neighbors"] union {w};
                    vertexData[ w ]["neighbors"]:=vertexData[ w ]["neighbors"] union {v};    
                    E := E union {{v,w}};
                    if not assigned(pathData[sites[1]]) then
                        pathData[sites[1]]:=table();
                        pathData[sites[1]]["vCell"]:={};
                    fi;
                    if not assigned(pathData[sites[2]]) then
                        pathData[sites[2]]:=table();
                        pathData[sites[2]]["vCell"]:={};
                    fi;
                    pathData[sites[1]]["vCell"] :=pathData[sites[1]]["vCell"] union {w,v};
                    pathData[sites[2]]["vCell"] :=pathData[sites[2]]["vCell"] union {w,v};
                
                fi; 
            fi;       
        od;
        mm:=mm minus {v}; 
        vertexData[ v ]["sites"]:=vSites[1];
    od;
    vertexData[ mm[1] ]["sites"]:=vertexData[ mm[1] ]["sites"][1];


         
    if minDigits <  ceil(-log10(infDist)) then  
        UI_LDigits := UI_LDigits  + ceil(-log10(infDist))-minDigits; #6 or more
        UI_HDigits := UI_Digits + ceil(-log10(infDist))-minDigits; 
        UI_Digits := UI_HDigits + ceil(-log10(infDist))-minDigits;          
        minDigits :=ceil(-log10(infDist));
        if  member('bug', _Env_GlobalOptions ) then
            print(cat("the estimated infimal distance has droped from ", 
            convert(infDist,string)," to ", convert(m,string)));
            print( "Working Digits are currently", UI_LDigits,UI_Digits,UI_HDigits);
        fi;            
    fi;       
end proc:

intersectionOfBLines:=proc(s1,s2,s3)
    local x1:=Re(s1), y1:=Im(s1), 
    x2:=Re(s2), y2:=Im(s2), 
    x3:=Re(s3),  y3:=Im(s3),
    m1,m2,m3,m4,r12,r13,r23,v,
    midpt12:=(x1+x2)/2+I*(y1+y2)/2,
    midpt23:=(x2+x3)/2+I*(y2+y3)/2,
    midpt13:=(x1+x3)/2+I*(y1+y3)/2;
    Digits := UI_HDigits;
    if member('bugV2', _Env_GlobalOptions ) then
        print("computing intersection for", s1,s2,s3);
        if toFloat(y1-y3) = 0 then
            r13:=I;
            
        else 
            r13:=toFloat((x1-x3)/(y3-y1));
        fi; 
    fi;    
    if toFloat(y1-y2) = 0 then
        r12:=I;
        m2:=toFloat(midpt12 + I*boundaryPoints["radius"]);   
    else 
        r12:=toFloat((x1-x2)/(y2-y1));
        m2:=toFloat(midpt12 + I*r12*boundaryPoints["radius"]+boundaryPoints["radius"]);
    fi;
    if toFloat(y2-y3) = 0 then
        r23:=I;
        m3:=toFloat(midpt23 + I*boundaryPoints["radius"]);  
    else 
        r23:=toFloat((x2-x3)/(y3-y2));
        m3:=toFloat(midpt23 + I*r23*boundaryPoints["radius"]+boundaryPoints["radius"]);
    fi;   
   
    #https://en.wikipedia.org/wiki/Line%E2%80%93line_intersection
    m1:=toFloat( midpt12 );
    m4:=toFloat( midpt23 );
    if member('bugV2', _Env_GlobalOptions ) then
        m5:=toFloat(midpt13 - ( I*r13*boundaryPoints["radius"]+boundaryPoints["radius"]));
        m6:=toFloat(midpt13 + ( I*r13*boundaryPoints["radius"]+boundaryPoints["radius"]));
        print("data");
        print("slope between s1,s2=",r12);  
        print("deltay=",y1-y2 );
         print("deltax=",x2-x1 );
        print("slope between s2,s3=",r23);
        print("deltay=",y2-y3 );
        print("deltax=",x3-x1 );
        print("two points on line 23");
        print(
            plots[display]({
                plot({complexLToPtL([m1,m2]),
                    complexLToPtL([m3,m4]),complexLToPtL([m5,m6])},color="green"),
                plots[pointplot]({complexToPt(midpt12),complexToPt(midpt23)},color="green"),    
                plots[pointplot](complexToPt(s1),color="Red"),
                plots[pointplot](complexToPt(s2),color="blue"),
                plots[pointplot](complexToPt(s3),color="black")
                }
                ,axes=none,title="s1=red,s2=blue,s3=black"));
    fi;
    v:=lineIntersection(m1,m2,m3,m4);
    
    if member('bugV3', _Env_GlobalOptions ) then
        m5:=toFloat(midpt13 - ( I*r13*boundaryPoints["radius"]+boundaryPoints["radius"]));
        m6:=toFloat(midpt13 + ( I*r13*boundaryPoints["radius"]+boundaryPoints["radius"]));
        print("data");
        print("slope between s1,s2=",r12);  
        print("deltay=",y1-y2 );
        print("deltax=",x2-x1 );
        print("slope between s2,s3=",r23);
        print("deltay=",y2-y3 );
        print("deltax=",x3-x1 );
        print("two points on line 23");
        print(
            plots[display]({
                plot({complexLToPtL([m1,m2]),
                    complexLToPtL([m3,m4]),complexLToPtL([m5,m6])},color="green"),
                plots[pointplot]({complexToPt(midpt12),complexToPt(midpt23)},color="green"),    
                plots[pointplot](complexToPt(s1),color="Red"),
                plots[pointplot](complexToPt(s2),color="blue"),
                plots[pointplot](complexToPt(s3),color="black"),
                plots[pointplot](complexToPt(v),color="red")
                }
                ,axes=none,title="s1=red,s2=blue,s3=black"));
    fi;
    return v;
end proc;

produceHurwitz := proc(); 
    local i,j, si,siteIndex, crtPt,startPt, s1OnRight,s1State,s1Dir,
    sitesInCommon,possDir,siWasOnRight,sitesOfLoc , nxtPt, 
    p1, p2, p3, crtNeighbors, mm,
    prevPt;

    Digits:=UI_LDigits;
    if member('bugH',_Env_GlobalOptions ) or member('bug',_Env_GlobalOptions ) then 
        print("~~~~~H");
       
        print(cat("Computing Hurwitz of ", convert(numOfProblemPts,
        string), " points."));
        print("all the sites");
        print(problemPoints);
        print(cat("Minimal distance between vertex,site collection is around ", convert( infDist,string)));
        print(cat("basepoint ", convert( basePoint,string)));
    fi;
    for siteIndex from 1 to  numOfProblemPts do 
        #for each site in order        
        crtPt := basePoint;   #start from canonnical vertex
        prevPt := NULL;
        si:=problemPoints[siteIndex];
        sitesOfLoc := vertexData[crtPt]["sites"];
        p1:=[];       
        siWasOnRight:=NULL;
        arrivedOnBoundary:=false;
        if member('bugH',_Env_GlobalOptions ) or  member('bug',_Env_GlobalOptions ) then 
            print("~~~~~ BEGINING Phase1 ~~~~~");
            print(cat(" PATH FOR ", convert(si=getSiteIndex(si),
            string)));
            print("sites of loc", sitesOfLoc);
        fi;
        while not member(si, sitesOfLoc) do 
            crtNeighbors := vertexData[crtPt]["neighbors"]
                 minus { prevPt,crtPt };
            if nops(crtNeighbors)>1 then
                crtNeighbors:=
                    sortArg(convert(crtNeighbors,list), crtPt); 
                s1Sites := [ seq(                         
                                ifelse( member( infinity,  vertexData[crtPt]["sites"] intersect vertexData[crtNeighbors[i]]["sites"]),
                                    [((vertexData[crtPt]["sites"] intersect vertexData[crtNeighbors[i]]["sites"]) minus{infinity})[1],infinity],
                                    sortArg(convert( vertexData[crtPt]["sites"] intersect vertexData[crtNeighbors[i]]["sites"],list),basePoint))
                            ,i=1..nops(crtNeighbors))];
                s2Sites := [ seq( (s1Sites[i])[2],i=1..nops(crtNeighbors))];   
                s1Sites := [ seq( (s1Sites[i])[1],i=1..nops(crtNeighbors))];     
                s1Dirs := [seq(normal(convert( s1Sites[i],rational) - convert(crtPt,rational)),i=1..nops(crtNeighbors))];
                s1Vectors:=[seq(normal(convert(crtNeighbors[i],rational) - convert(crtPt,rational)),i=1..nops(crtNeighbors))];  
                toRights:=[seq( wIsToRightOfz(s1Vectors[i],s1Dirs[i]),i=1..nops(crtNeighbors))];   
                s1States:=[seq( 
                    ifelse( siteIndex > getSiteIndex(s1Sites[i]),"OLD","NEW"),
                    i=1..nops(crtNeighbors))];  
                s2States:=[seq( 
                    ifelse( siteIndex > getSiteIndex(s2Sites[i]),"OLD","NEW"),
                    i=1..nops(crtNeighbors))];                                
                
                if member('bugH',_Env_GlobalOptions ) then
                    print("~~~~~~~~~~"); 
                    print(cat("Site: ",convert(si, string)));
                    print(cat("Crt location: ", convert(crtPt, string)));
                    print(cat("Arrived on boundary?", convert(arrivedOnBoundary,string)));
                    print(cat("On boundary? ", convert(vertexData[crtPt]["onBoundary"],string)));
                    print(cat("nearby sites: ",convert(sitesOfLoc,string)));    
                    print(cat("current Niehbors: ",convert(crtNeighbors,string)));
                     print(cat("youngest Sites: ", convert(s1Sites, string)));
                    print(cat("age of Sites: ", convert(getSiteIndex~(s1Sites), string))); 
                    print(cat("s1states: ", convert(s1States, string)));
                     print(cat("youngest to right? ",convert(toRights,string))); 
                     
                     print(cat("oldest Sites: ", convert(s2Sites, string)));
                    print(cat("age of Sites: ", convert(getSiteIndex~(s2Sites), string)));                     
                    print(cat("s2 states ",convert(s2States,string)));
                   
                    if member('bugHVisual',_Env_GlobalOptions ) then
                        print(plots[display]( { plots[pointplot](
                            {complexToPt(si)},color="red",symbol=cross),
                            plots[pointplot](
                            {complexToPt(s1Sites[1])},color="green",symbol=cross),
                            plots[pointplot](
                            {complexToPt(s2Sites[1])},color="orange",symbol=cross),
                        ifelse(p1<>[],plot(complexLToPtL~(p1), color="black"),NULL),
                        plots[pointplot](complexToPt(crtNeighbors[1]),color="green"),
                        plots[pointplot](complexToPt(crtNeighbors[-1]),color="orange"),
                        ifelse(nops(complexToPt)>2, 
                            plots[pointplot](complexToPt~(crtNeighbors[2..-2]),color="blue"),NULL),
                        plots[pointplot]({complexToPt(crtPt)}, color="green"),
                        ifelse(siteIndex < numOfProblemPts-1,
                            plots[pointplot](complexToPt~( { seq(problemPoints[i], i=siteIndex+1..numOfProblemPts-1) }),color="black",symbol="cross"),NULL),
                        ifelse( siteIndex >1, 
                            plots[pointplot](complexToPt~( { seq(problemPoints[i], i=1..siteIndex-1) }),color="yellow"),NULL)},
                        title="Current state: target site is red, yellow are the other sites already attained, orange are the sites to be attained"));
                    fi;
                          
                fi;
                
                if s2Sites[1]=infinity and s1States[1]= "NEW" then
                    II:=1;
                elif s2Sites[-1]=infinity and s1States[-1]= "OLD" then
                    II:=-1;
                else
                    II:=1;
                    while s2States[II] = s1States[II]  do
                        II:=II+1;
                    od;
                    if not toRights[II] then
                        II:=II+1;
                    fi;
                fi;
                #if (toRights[II] and s1States[II]="OLD" and s2States[II]="NEW") or (not toRights[II] and s1States[II]="NEW" and s2States[II]="OLD" ) then
                #    II:=II;
               # elif (toRights[II+1] and s1States[II]="OLD" and s2States[II]="NEW") or (not toRights[II] and s1States[II]="NEW" and s2States[II]="OLD" ) then
                #    II:=II;
                #    II:=II-1;
                #fi;    
            else
                if member('bugH',_Env_GlobalOptions ) then
                    print("~~~~~~~~~~"); 
                    print("only one nbr");                          
                fi;
                II:=1;
            fi;       
            arrivedOnBoundary:=vertexData[crtNeighbors[II]]["onBoundary"] and not vertexData[crtPt]["onBoundary"];
            sitesOfLoc := vertexData[ crtNeighbors[II]]["sites"];  
            p1 := [ op(p1) , [ crtPt , crtNeighbors[II]]];
            if not assigned(edgeData[{crtPt,crtNeighbors[II]}]) then
                edgeData[{crtPt,crtNeighbors[II]}]:=table();
                edgeData[{crtPt,crtNeighbors[II]}]["normDir"] := [crtPt , crtNeighbors[II]]; 
                edgeData[{crtPt,crtNeighbors[II]}]["gam"] :=normal(convert((crtNeighbors[II]-crtPt)/2,rational));
                edgeData[{crtPt,crtNeighbors[II]}]["path"] :=crtPt*(1-s)/2+crtNeighbors[II]*(1+s)/2;
            fi; 
            prevPt:=crtPt;    
            crtPt := crtNeighbors[II];
            if not member(si, sitesOfLoc) then                                       
                siWasOnRight := s1OnRight;           
            fi;    
        od;
        prevPt:=NULL;
        siWasOnRight:=NULL;
        s1OnRight:=NULL;        

        p2:=[];
        startPt:=crtPt;
        
        if p1<>[] and (member('bugH',_Env_GlobalOptions ) or member('bug',_Env_GlobalOptions )) then
            print("~~~~~Phase1 construction complete~~~~~");
            print(plots[display]( 
                { plots[pointplot](
                    {complexToPt(si)},color="red",symbol=diamond),
                plot(complexLToPtL~(p1), color="black"),
                ifelse(siteIndex < numOfProblemPts-1,
                    plots[pointplot](complexToPt~( { seq(problemPoints[i], i=siteIndex+1..numOfProblemPts-1) }),color="black",symbol=diamond),NULL),
                ifelse( siteIndex >1, 
                    plots[pointplot](complexToPt~( { seq(problemPoints[i], i=1..siteIndex-1) }),color="yellow"),NULL)
                    },
                title="Current state: target site is red, yellow are the other sites already attained, orange are the sites to be attained"));    
            print("~~~~~Phase2 construction ~~~~~");
        fi;
             
        while crtPt <> startPt or prevPt=NULL do
            crtNeighbors := vertexData[crtPt]["neighbors"] intersect
                pathData[si]["vCell"] minus { prevPt };            
            crtNeighbors:=op~(constructEdge(crtNeighbors,crtPt) intersect E) minus {crtPt};
            if member('bugH',_Env_GlobalOptions ) then
                print(plots[display]( { plots[pointplot](
                    {complexToPt(si)},color="red",symbol=diamond),
                ifelse(p1<>[],plot(complexLToPtL~(p1), color="black"),NULL),
                ifelse(p2<>[],plot(complexLToPtL~(p2), color="black"),NULL),
                plots[pointplot](complexToPt(crtNeighbors[1]),color="green"),
                plots[pointplot](complexToPt(crtNeighbors[-1]),color="orange"),
                ifelse(nops(complexToPt)>2, 
                    plots[pointplot](complexToPt~(crtNeighbors[2..-2]),color="blue"),NULL),
                plots[pointplot]({complexToPt(crtPt)}, color="green"),
                ifelse(siteIndex < numOfProblemPts-1,
                    plots[pointplot](complexToPt~( { seq(problemPoints[i], i=siteIndex+1..numOfProblemPts-1) }),color="black",symbol=diamond),NULL),
                ifelse( siteIndex >1, 
                    plots[pointplot](complexToPt~( { seq(problemPoints[i], i=1..siteIndex-1) }),color="yellow", symbol=diamond),NULL)},
                title="Current state: target site is red, yellow are the other sites already attained, orange are the sites to be attained"));    
            fi;
            if nops(crtNeighbors)>1 then                
                crtNeighbors:=sortArg(convert(crtNeighbors,list), crtPt); 
                s1Sites := [ seq(                         
                                ifelse( 
                                    member( infinity,  vertexData[crtPt]["sites"] intersect vertexData[crtNeighbors[i]]["sites"]),
                                    [op((vertexData[crtPt]["sites"] intersect vertexData[crtNeighbors[i]]["sites"]) minus{infinity}),infinity],
                                    sortArg(convert( vertexData[crtPt]["sites"] intersect vertexData[crtNeighbors[i]]["sites"],list),basePoint))
                            ,i=1..nops(crtNeighbors))];
                s2Sites := [ seq( (s1Sites[i])[2],i=1..nops(crtNeighbors))];   
                s1Sites := [ seq( (s1Sites[i])[1],i=1..nops(crtNeighbors))];     
                s1Dirs := [seq(normal(convert( s1Sites[i],rational) - convert(crtPt,rational)),i=1..nops(crtNeighbors))];
                s1Vectors:=[seq(normal(convert(crtNeighbors[i],rational) - convert(crtPt,rational)),i=1..nops(crtNeighbors))];  
                toRights:=[seq( wIsToRightOfz(s1Vectors[i],s1Dirs[i]),i=1..nops(crtNeighbors))];   
                
                II:=1;

                while  (getSiteIndex(s1Sites[II])< siteIndex and getSiteIndex(s2Sites[II])< siteIndex) or
                    (s1Sites[II]=si and toRights[II]) or
                     (s2Sites[II]=si and not toRights[II]) do
                        II:=II+1;                        
                od;
                if member('bugH',_Env_GlobalOptions ) then
                    print("~~~~~Nbr identification~~~~~"); 
                    print(cat("nbrs: ",convert(crtNeighbors ,string))); 
                    print(cat("finite sites: ", convert( s1Sites, string)));
                    print(cat("s1 is to the right of path? ",convert(toRights,string)));                                                                          
                    print(cat("index ", convert(II, string)))
                fi;
                crtNeighbors:={crtNeighbors[II]};   
            fi;             
            
                      
            p2 := [ op(p2) , [ crtPt , crtNeighbors[1]]];
            if not assigned(edgeData[{crtPt,crtNeighbors[1]}]) then
                edgeData[{crtPt,crtNeighbors[1]}]:=table();
                edgeData[{crtPt,crtNeighbors[1]}]["normDir"] := [crtPt , crtNeighbors[1]]; 
	            edgeData[{crtPt,crtNeighbors[1]}]["gam"] :=normal(convert((crtNeighbors[1]-crtPt)/2,rational));
                edgeData[{crtPt,crtNeighbors[1]}]["path"] :=crtPt*(1-s)/2+crtNeighbors[1]*(1+s)/2;
            fi; 
            prevPt:=crtPt;           
            crtPt:=crtNeighbors[1];          
        od;  
        if member('bugH',_Env_GlobalOptions ) then
            print("~~~~~Phase2 construction complete~~~~~"); 
            print(plots[display]( { plots[pointplot](
                {complexToPt(si)},color="red",symbol=cross),
                ifelse(p1<>[],plot(complexLToPtL~(p1), color="black"),NULL),
                ifelse(p2<>[],plot(complexLToPtL~(p2), color="black"),NULL),
                ifelse(siteIndex < numOfProblemPts-1,
                plots[pointplot](complexToPt~( { seq(problemPoints[i], i=siteIndex+1..numOfProblemPts-1) }),color="ORANGE"),NULL),
                ifelse( siteIndex >1, 
                plots[pointplot](complexToPt~( { seq(problemPoints[i], i=1..siteIndex-1) }),color="yellow"),NULL)},
                title="Current state: target site is red, yellow are the other sites already attained, orange are the sites to be attained"));
        fi;    
        if nops(p1)>0 then           
            p3 := reverseList(map(reverseList,p1))
        else
            p3:=[];
        fi;        
        pathData[si]["path"] :=  [op(p1),op(p2), op(p3)]; 
        if member('visual',_Env_GlobalOptions )
        or member('bug',_Env_GlobalOptions ) or member('showpaths',_Env_GlobalOptions ) then 
            updatePlotter( # [{ifelse( si < numOfProblemPts, [problemPoints[si]], boundaryPoints["boundary"]),  {color=colors[(si mod nops(colors))+1],symbol=asterisk}}], 
                    [{}],[{complexLToPtL~(pathData[si]["path"]), {color="red"}}
                ]);  
        fi;        
        if member('bugH1',_Env_GlobalOptions )  then 
            print(cat("path has been completed for: ",
            convert(si,string))); 
            print("with the first portion of path:");
            print(p1); 
            print("with the second portion of path:");
            print(p2);
            print("with the third portion of path:");
            print(p3);
           
        fi; 
    od;
end proc;   

generateMonodromy := proc(x::name, y::name)
    local
	i, j, k,l, pt1,pt2,fibre, yi, yj, yk,sigma_inf,sigma,
	siteToMonoData:=Array([]),mm,val,
	si,fibre1,fibre2,edgesBlack:={}, edgesGreen:={};  
    if member('bug',_Env_GlobalOptions )  then 
        print("now obtaining a geometric monodromy"); 
        print("edges:");
        #print(E);                         
    fi;
    
    Digits:=UI_LDigits;  
    for e in E do                 
        aContEdge(e,x,y); 
        if member('visual',_Env_GlobalOptions) then 
            if  convert(edgeData[e]["perm"], disjcyc) <> convert(Perm[[]], disjcyc) then           
                edgesGreen:=edgesGreen union {complexLToPtL([op(e)])};
            else
                edgesBlack:=edgesBlack union {complexLToPtL([op(e)])}; 
            fi;            
        fi;
    od; 
    if member('visual',_Env_GlobalOptions) then  
        print(plots[display]({plot(edgesBlack,color="black"),plot(edgesGreen,color="green")},axes=none,title="Edges that permute fibres non trivially are green"));                
    fi; 
            
    singPoints:=Array([]);   
    branchPoints:=Array([]);   
    HurwitzFormula:=0;
    for siteIndex from 1 to numOfProblemPts do  #for each path
        si:=problemPoints[siteIndex];
        if member('bugPerm',_Env_GlobalOptions) then
            print(cat("site",convert(siteIndex ,string)));
            print(cat("has ", convert(nops(pathData[si]["path"]),string), " edges"));  
            print(pathData[si]["path"]);                      
        fi;
        pathData[si]["perm"]:= Perm([]);   
        for e in pathData[si]["path"] do
            if isNormalForm(e) then
                sigma_e := edgeData[{op(e)} ]["perm"]; 
                if member('bugPerm',_Env_GlobalOptions) then
                   
                    print(e," is normal form");
                    print("AC");
                    print(edgeData[{op(e)}]["AC"]);
                    print("First fib should be ordered by lex");
                    print(edgeData[{op(e)}]["AC"][1][2]);
                    print("Final fib");
                    print(edgeData[{op(e)}]["AC"][-1][2]);
                    print("is a permuation of lex order");
                    print(sortLex(edgeData[{op(e)}]["AC"][-1][2]));
                    print(getOrderedRoots(x,y,F[1],e[2]));
                    print("perm", sigma_e);
                fi; 
            else 
                mm:=reverseEdgeAC(e);
                sigma_e:=mm[1];
                if member('bugPerm',_Env_GlobalOptions) then
                    print(e," is not normal form");
                    print("AC");
                    print(mm[2]);
                    print("First fib should be ordered by lex");
                    print(mm[2][1][2]);
                    print("Final fib");
                    print(mm[2][-1][2]);
                    print("is a permuation of lex order");
                    print(sortLex(mm[2][-1][2]));
                    print(getOrderedRoots(x,y,F[1],e[1]));
                    print("perm", sigma_e);
                fi;                
            fi;
            #travel f to g then produce is fg
            pathData[ si]["perm"]:= 
                pProd(pathData[si]["perm"],sigma_e); 
            if member('bugPerm',_Env_GlobalOptions) then
                print("the permutation is updated to");
                print(pathData[ si ]["perm"]);
            fi;                
        od;    

        if pathData[si]["perm"] <> Perm([]) then  
            permRep:=convert(pathData[si]["perm"] , 'disjcyc');          
        	pathData[si]["NumOfPlaces"]:=nops(permRep);
            pathData[si]["Branching"]:=[nops~(permRep),permRep]; #perhaps can be used to determine puisuex
            
            Digits:=UI_ODigits;
            ArrayTools[Append](siteToMonoData,[
                ifelse(si <>infinity, 
                    scrub(si), 
                   si),
                    permRep 
                    ] );       
            ArrayTools[Append](branchPoints,si );            
        else 
            pathData[si]["NumOfPlaces"]:=numSheets;      
			ArrayTools[Append](singPoints,si );        
        fi;  
    od;
    sigma_e:=Perm([]);
    for i from 1 to numOfProblemPts do

        sigma_e:=pProd(sigma_e,pathData[problemPoints[i]]["perm"] );        
    od;
    
    #if HurwitzFormula/2-numSheets+1<>g then
    #    print(HurwitzFormula/2-numSheets+1);
    #    error "check Monodromy, Hurwitz not satisfied"
    if sigma_e <> Perm([]) then
        print(sigma_e);
        clearCanvas();  
        error "the product of the monodromy generators is not the identity";
    fi;
    #get formula for another check
	#check_cstruct( pathData[numOfProblemPts]["perm"] ,algcurves[puiseux](F[1],x=infinity,y,`cycle structure`));
   
    siteToMonoData;       
end proc:
check_cstruct:=proc(perm,val)
	local a;
	a:=map(nops,convert(perm,'disjcyc'));   
		if sort([1$(convert(val,`+`)-convert(a,`+`)),op(a)])<>sort(val) then
			clearCanvas();  
            error "Found wrong cycle structure"
		fi
end:



removeVerts := proc(pt1)
    local pt2,nbrs;
        nbrs:= vertexData[pt1]["neighbors"];
        for pt2 in nbrs do
            
            E:=E minus {{pt2,pt1}};
            print("SSFE");
            print({pt2,pt1});
            vertexData[pt2]["neighbors"]:=vertexData[pt2]["neighbors"] minus{pt1};
        od;
        V:=V minus {pt1};
        for site in vertexData[pt1]["sites"] do
            pathData[site]["vCell"]:= pathData[site]["vCell"] minus {pt1};
        od;
        unassign(vertexData[pt1]);
        gc();

end proc: