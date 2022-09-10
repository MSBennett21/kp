# PeriodMatrix
# DESCRIPTION
# working idea: compute mono to get hurwitz, get linear combos from homology, get the associated hurwitz eles and parametrize them. Then use intgrate comand against the differentials and Y as is.

periodMatrix := proc(curve, x::name, y::name)
	local pm;
	_Env_GlobalOptions := ifelse( assigned(_Env_GlobalOptions),
        [ op(_Env_GlobalOptions), args[4..-1], 'int', 'secondary' ],
		[args[4..-1], 'int', 'secondary']); 

	if member('bug', _Env_GlobalOptions) or member('bugI0', _Env_GlobalOptions) then
		print("Initalization phase of integration");
		print("global options are");
		print(_Env_GlobalOptions);
		print("calling homology");
	fi;
	UI_ODigits:=Digits;
	Digits:=UI_ODigits;
	Homology(curve,x,y);

	if member('bug', _Env_GlobalOptions) or member('bugI0', _Env_GlobalOptions) then
		print("Initalization phase of integration complete");
		print("The three sets of digits are now:");
		print(UI_LDigits,UI_Digits,UI_HDigits);
		print("with Digits being ",Digits);
		print("and user specified ", UI_ODigits);
	fi;
	
	pm := periodMatrixInner(x,y); #it would be good to pass over monodromy data for faster speeds.	
	clearCanvas();
	pm;
end:

periodMatrixInner := proc(x::name, y::name)
	local pm, A, B,i,j,m;
	pm,A,B := Pmatrix(x,y,diffBasis); #it would be good to pass over monodromy data for faster speeds.
		
	if member('R',_Env_GlobalOptions) then
		if member('bug', _Env_GlobalOptions) or  member('bugI0', _Env_GlobalOptions) then
			print("period matrix before processing");
			print(pm);
			print("A:");
			print(A);
			print("B:");
			print(B);
		fi;	
		pm := linsolve(A,B,method=none,
			inplace=false, outputoptions=[]); 
		if g>1 then
			m:=max( seq(seq( abs(pm[i,j]-pm[j,i])/2, i=j+1..g), j=1..g-1) ); 
			i:=cat("About: ",convert(evalf(-log10(m),3),string),
				" accurate Digits");
			if m > 10^(3-UI_ODigits) then # there are entries with digits diffing more than digits-3
				WARNING( i ) # Shouldn't happen.
			else
				userinfo(2,'algcurves',i) #look up user info
			fi;
			if member('bug', _Env_GlobalOptions) or  member('bugI0', _Env_GlobalOptions) then
				print("period matrix is");
				print(pm);
				print("Inf norm of the skew symetric part: ", m);
			fi;	
		fi;
		# Test if Imaginary part is positive definite.
		m:=Im(pm);
        #print(m);
		if not LinearAlgebra[IsDefinite](m,query='positive_definite') then
			clearCanvas();  
			error "Imaginary part not positive definite" 
		fi;
		
		pm	
	else
		if member('bug', _Env_GlobalOptions) or member('bugI0', _Env_GlobalOptions) then
			print("period matrix after processing");
			print(pm);
		fi;	
		pm 
	fi;
end:
#use crtsheet nxt sheet when on negs, term and start when on lowlevel
Pmatrix:=proc(x::name,y::name,f::Vector)
	local vecVal,crtSheet,A,pm, B,i,j,nxtSheet,
	sI,cj,mm,cycleData; 
	
	for i from 1 to 2*g do				
		cycleData:=op(homBasis[i]);
		if member('bug', _Env_GlobalOptions) or member('bugI0', _Env_GlobalOptions)  then
			print(cat("CONSIDERING HOMOLOGY ELEMENT ",
			ifelse(i<=g, 
				cat("a",convert(i,string)),
				cat("b",convert(i,string)))
			));
			
			#print("what we recieve from hom",homBasis[i]);
			print("relevant data");
			print(cycleData);	
		fi;	
		vecVal:=prod(0,f);
		for j from 1 to nops(cycleData)/2 do	
			cj:=cycleData[2*j][2]; #cycle cj times
			sI:=getSiteIndex(cycleData[2*j][1]); #around site I				
			
			if member('bugI', _Env_GlobalOptions) then
				print("~~I0: HOM element Data");
				print(cat("integrating around ", convert(sI=cycleData[2*j][1],string))); 
				print(cat(convert(cj,string), "times"));
			fi;	
			#FACT if it worked we will have traveled from crt sheet to next sheet.
			crtSheet := cycleData[2*j-1];
			nxtSheet := ifelse(2*j+1 < nops(cycleData), cycleData[2*j+1],cycleData[1]); 
			if cj < 0 then
				# if cj is negative then we are actually starting at next sheet and ending at crt sheet
				if member('bugI', _Env_GlobalOptions) then
					print(cat("TERMINAL SHEET I0: ", convert( crtSheet,string)));
					print(cat("START SHEET I0: ", convert(nxtSheet,string)));
					print("cj is negative and so we are integrating backwards") ;
					print("I0 end  ~~");
				fi;	
				mm := intBpNTimes(x,y,abs(cj),getSite(sI),nxtSheet,f);	
				if mm[2]<> crtSheet then
					clearCanvas();
					error "sheets not aligned:BP11";
				fi;
			else
				if member('bugI', _Env_GlobalOptions) then
					print(cat("TERMINAL SHEET I0: ", convert( nxtSheet,string)));
					print(cat("START SHEET I0: ", convert(crtSheet,string)));
					print("I0 end  ~~");
				fi;
				mm := intBpNTimes(x,y,cj,getSite(sI),crtSheet,f);	
				if mm[2]<> nxtSheet then
					clearCanvas();
					error "sheets not aligned:BP12";
				fi;
			fi;
					
			if member('bugI', _Env_GlobalOptions) then
				print("~~I0: values ");
				print("sign of Integral: ", cj/abs(cj));
				print("result");
				print(mm);
				print("I0 end  ~~");
			fi;	
			vecVal:=vecVal+prod(cj/abs(cj),mm[1]);
		od;	
		if i<=g	then 
			A:=ifelse(i=1, vecVal,<A|vecVal>);
		else
			B:=ifelse(i=1+g, vecVal,<B|vecVal>);
		fi;	
		if member('bug', _Env_GlobalOptions) or member('bugI0', _Env_GlobalOptions) then
			print("~~Matrix Construction : current state ");
			print("A: ");print(A);
			print("B:");print(B);
		fi;	
	od;
	pm:=<A|B>;
	return pm, A, B;
end:
intBpNTimes:= proc(x::name,y::name, cj::integer,site, stSheet::integer,f::Vector)
	local i,val, vecVal:=prod(0,f),mm, sigma:=pathData[site]["perm"],
	startSheet:=stSheet,termSheet:=sigma[startSheet];
	if cj<=0 then error "integrator integrates in one direction only"; fi;
	if member('bugI1', _Env_GlobalOptions) then
		print("~~I1: BP INTEGRATION LVL 1");
		print("START SHEET I1: ", stSheet);
		print("END SHEET (THEO) I1: ", termSheet);
	fi;
	for i from 1 to cj do #for each of these
		if member('bugI1', _Env_GlobalOptions) then
			print("~~I1: iterations");
			print("NUM OF TIMES AROUND SITE SO FAR: ", i);
			print("VALS: ", vecVal);
			print("I1~~~");
		fi; 	
		mm:=intAroundBp( x, y, site,startSheet, f );
		if member('bugI1', _Env_GlobalOptions) then
			print("~~I1: UPDATE");
			print("CRT VALS:", mm);
		fi; 
		startSheet:=mm[2];
		vecVal:=vecVal+mm[1];
		if startSheet<> termSheet then
			clearCanvas();
			error cat("The permutation around site did not match the outcome of the integraton:I1");
		fi;
		termSheet:=sigma[startSheet];
	od;
	
	return [vecVal,startSheet];
end proc;

intAroundBp:= proc(x::name,y::name,site,stSheet::integer,f::Vector)
	local vecVal,edge,crtSheetEdge,sigma_e,
	nxtSheetEdge,e,eIsNormalForm;
	#print(sI,stSheet);
	if member('bugI2', _Env_GlobalOptions) then
		print("~~I2: BP INTEGRATION LVL 2");
		print("START SHEET I2: ", stSheet);
	fi;	
	if assigned(pathData[site]["int"][stSheet]) then	
		if member('bugI2', _Env_GlobalOptions) then
			print(cat("This integral has already been stored I2"));
		fi;	
		return pathData[site]["int"][stSheet]; #if we already have the bp for this sheet just return it
	else
		if member('bugI2', _Env_GlobalOptions)then
			print(cat("This integral must be formed I2."));
			print("PATH OF INTEGRATION");
			print(pathData[site]["path"]);
		fi;	
		vecVal:=prod(0,f);
		crtSheetEdge:=stSheet;
		for edge in pathData[site]["path"] do		#otherwise we must integrate edge by edge starting at the ref sheet
			sigma_e:=edgeData[{op(edge)}]["perm"];
			
			if member('bugI2', _Env_GlobalOptions) then	
				print("~~I2: EDGE INTEGRATION LVL 2");
				print(cat("Integrating ",
				convert(edge,string)));	
				print("NORMAL FORM ",edgeData[{op(edge)}]["normDir"]);			
			fi;	
			if isNormalForm(edge) then
				e:=edge;
				nxtSheetEdge:=sigma_e[crtSheetEdge];
				eIsNormalForm:=1;
				if member('bugI2', _Env_GlobalOptions) then
					print("START SHEET I2", crtSheetEdge);
					print("Terminal SHEET (THEO) I2: ", nxtSheetEdge);
					print("I2~~~");
				fi;	
				mm:=intEdge(x,y,e,crtSheetEdge,f);
				if nxtSheetEdge<>mm[2] then
					clearCanvas();
					error 
					cat("Edge did not integrate from",convert(crtSheetEdge,string),
					" to ", convert( nxtSheetEdge, string)," I21");
				fi;
				crtSheetEdge:=nxtSheetEdge;
			else
				e:=[edge[2],edge[1]];
				nxtSheetEdge:=pInv(sigma_e)[crtSheetEdge];
				eIsNormalForm:=-1;
				if member('bugI2', _Env_GlobalOptions) then
					print("START SHEET I2: ", nxtSheetEdge);
					print("Terminal SHEET (THEO) I2: ", crtSheetEdge);
					print("DIRECTION IS NEGATIVE");
					print("I2~~~");
				fi;
				mm:=intEdge(x,y,e,nxtSheetEdge,f);
				if crtSheetEdge<>mm[2] then
					print("!!!!");
					print(mm);
					error 
					cat("Edge did not integrate from",convert(nxtSheetEdge,string),
					" to ", convert( crtSheetEdge, string), " I22");
				fi;
				crtSheetEdge:=nxtSheetEdge;
			fi;
			
			vecVal := vecVal+prod(eIsNormalForm,mm[1]);

			if member('bugI2', _Env_GlobalOptions) then
				print(cat("Integration of edge complete. Updated value: "));	
				print(vecVal);
				print("I2~~~");
			fi;
		od;	
		if member('bugI2', _Env_GlobalOptions) then
			print("~~~~~BP2");
			print(cat("Integration of bp complete. Updated value: "));	
			print(vecVal);
			print("arrived to sheet ", crtSheetEdge);
			print("I2~~~");
		fi;						
		pathData[site]["int"][stSheet]:=[vecVal,crtSheetEdge];
		return pathData[site]["int"][stSheet];
	fi;
end proc;
intEdge:= proc(x::name,y::name,edge,stSheet::integer,f::Vector)
	local mlvl,mm1,mm2,mm3,
	AC:= edgeData[{op(edge)}]["AC"],
	sigma:=edgeData[{op(edge)}]["perm"],		
	termSheet :=sigma[stSheet],
	gam:=edgeData[{op(edge)}]["gam"],
	tol:=1, Lvl:=BLvl-1,vecVal:=prod(0,f),
	path:=edgeData[{op(edge)}]["path"],
	pathPos:=(edge[1]-edge[2])/2*s+edge[2],
	pathNeg:=-(edge[1]-edge[2])/2*s+edge[1],
	intAtLvl:= proc(e,sheet::integer,Lvl::integer)
		local eps:=1, II:=edgeData[{op(e)}]["zeroIndex"],
		k:=0,kk,xPos,yPos,xNeg,yNeg,xPivPos,xPivNeg,sPivPos,sPivNeg,
		mm,J:=2^(-Lvl), maxVal:=-1;
		#at this lvl, <-inner function
		#while each summand has magnitude greater than tolerance1
		Digits:=3*UI_Digits;
		if GLvl<Lvl then
			GLvl:=Lvl;
		fi;
		if member('bugI4', _Env_GlobalOptions) then
			print("~~I4 INTGRATION AT LVL");
			print("Global Lvl, ", GLvl);
			print(cat("Edge step h is: ", convert(J,string))); 	
		fi;	
		if Lvl=BLvl then
			wt:=intData["w"][0];	
			xPos:=normal(subs(s=(1-intData["s"][0]), path));
			yPos:=newtonY( subs(x=xPos,F[1]),y,
				AC[II][2][sheet]);
			mm:=evalf(subs(x=xPos,y=yPos,prod(wt,f)));	
		else
			mm:=prod(0,f);
		fi;
		while eps > 1/10^(3*UI_ODigits) do
			if member('bugI4', _Env_GlobalOptions) then
				print("~~I4 ITERATIONS"); 
				print("eps= ", eps);
				print("mm= ", mm);
			fi;
			k := k+1;
			kk := ifelse( Lvl=BLvl, k, 2*k-1);
			tVal := kk*J;
			if member('bugI4Lvl', _Env_GlobalOptions) then
				print("~~ I4 Lvl iteration 1");
				print("index =",kk);
				print(kk*J);
			fi;		
			if assigned(intData["s"][normal(tVal)]) then
				sVal := intData["s"][normal(tVal)];
				wt := intData["w"][normal(tVal)];
				if member('bugI4Lvl1',  _Env_GlobalOptions) then
					print("~~ I4 Lvl iteration 2");			
					print("sVal is already stored, lvl<=Glevel ");
				fi;			
			else
				wt:=evalf(subs(t=tVal,intData["Gt"]));
				sVal:=evalf(subs(t=tVal,intData["abs"]));
				intData["s"][normal(tVal)] := sVal;
				intData["w"][normal(tVal)] := wt;
				GN:=GN+1;
				if member('bugI4Lvl1',  _Env_GlobalOptions) then
					print("~~ I4 Lvl iteration 2");					
					print("sVal is not already stored, lvl=Glevel ");
					print("total number of absicussis",GN);
					print(1-intData["s"][normal(tVal)]);
					print(Lvl);
					print("wt");
					print(intData["w"][normal(tVal)]);
				fi;	
			fi;
			if k=1 then
				sPivPos:=convert(AC[II][1],rational);
				xPivPos:=normal(convert(subs(s=sPivPos,path),rational));
				yPivPos:=AC[II][2][sheet];
				sPivNeg:=convert(AC[II][1],rational);
				xPivNeg:=normal(convert(subs(s=sPivNeg,path),rational));
				yPivNeg:=AC[II][2][sheet];
			fi;
			if member('bugI4Lvl2',  _Env_GlobalOptions) then
				print("~~ I4 Lvl iteration 3");					
				print("pos Piv", xPivPos);
				print("pos yPiv", yPivPos);
				print("neg Piv", xPivNeg);
				print("neg yPiv", yPivNeg);
			fi;				
			xPos:=evalf(subs(s=sVal,pathPos));	
			xNeg:=evalf(subs(s=sVal,pathNeg));		
			der:=evalf(subs(x=xPivPos,evalf(subs(y=yPivPos, (xPos-x)*F["Yx"]+y))));						
			yPos:=newtonY(subs(x=xPos,F[1]),y,der);
			der:=evalf(subs(x=xPivNeg,evalf(subs(y=yPivNeg, (xNeg-x)*F["Yx"]+y))));				
			yNeg:=newtonY(subs(x=xNeg,F[1]),y,der);
			xPivPos:=xPos;
			xPivNeg:=xNeg;
			yPivPos:=yPos;
			yPivNeg:=yNeg;
			mPos:=evalf(
				subs(x=xPos,
						y=yPos,prod(wt,f)));
			mNeg:=evalf(
				subs(x=xNeg,
				y=yNeg,prod(wt,f)));
			if member('bugI5', _Env_GlobalOptions) then		
				print("Pos value");
				print(mPos);
				print("Neg value");
				print(mNeg);			
			fi;	
			eps:=max(norm(mPos,infinity),norm(mNeg,infinity));	
			maxVal:=max(eps,maxVal);
			mm:=mm+mPos+mNeg;		
		od;

		eps:=extractFib(yPos,
			sortLex(AC[-1][2]));
		return [mm,eps[2],maxVal];			
	end proc;
		

	if member('bugI3',_Env_GlobalOptions) then
		print("~~I3: EDGE INTEGRATION LVL2");
		print("START SHEET ", stSheet);
		print("PATH,", path);
		print("derivative, ", gam);
		print("THEORETICAL TERM, ", termSheet);
	fi;
	if not assigned(edgeData[{op(edge)}]["int"][stSheet]) then	
		if member('bugI3',_Env_GlobalOptions) then
			print("EDGE needs integration");
		fi;
		m3:=-1;crtMax :=-1;
		while tol > 1/10^(UI_ODigits) and Lvl<GlobalLvl do #while diff between prevVal-crtVal is large
			 #increase lvl of edge
			Lvl:=Lvl+1;
			if member('bugI3', _Env_GlobalOptions) then	
				print("~~I3: EDGE INTEGRATION itteration LVL2");		
				print(" CRT LVL: ",Lvl);
			fi;	
			mlvl :=intAtLvl(edge,stSheet,Lvl);
			crtMax :=max(crtMax,mlvl[3]);
			print(mlvl);
			vecVal := vecVal + mlvl[1];
			if Lvl=BLvl then
				m1 := prod(gam/2^(Lvl), vecVal);
				d1:=1;
				d2:=1;
				d3:=1;
				tol := 1;
			elif Lvl=BLvl+1 then
				m2:= prod(gam/2^(Lvl), vecVal);
				
			else
			
				m1:=ifelse(m3=-1, m1, m2);
				m2:=ifelse(m3=-1, m2, m3);
				m3:= prod(gam/2^(Lvl), vecVal);
				print(m1,m2,m3);
				d1:=log10(norm(m3-m2,infinity));
				d2:=log10(norm(m3-m1,infinity));
				d3:=log10(10^(-UI_Digits)*crtMax);
				tol:=10^(max(d1^2/d2,2*d1,d3));
			fi;
			
			if member('bugI33', _Env_GlobalOptions) then			
				print(" CRT tol: ",tol);
				print(" CRT d1 sn-sn-1: ",d1);
				print(" CRT d2 sn-sn-2: ",d2);
				print(" CRT d3 largest*eps: ",d3);

			fi;	
		od;
		edgeData[{op(edge)}]["int"][stSheet]:=[prod(gam/2^Lvl,vecVal),mlvl[2]];	
	fi;
	unassign(intAtLvl);
	return
		edgeData[{op(edge)}]["int"][stSheet];
end proc;


intProcessing := proc(x,y)

    diffBasis := algcurves[differentials]( F[1], x, y, 'skip_dx' ); 
	diffBasis := convert([seq(convert(diffBasis[i]/F["y"], horner,[x,y]),i=1..g)],Vector);
	if member('int', _Env_GlobalOptions) then
        print("processing differentials vecctor");
		print(diffBasis);
    fi;  
end proc;

#the integral on a line of a function f:C->C

	