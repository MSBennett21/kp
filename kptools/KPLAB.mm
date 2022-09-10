

RealKPSolLab := proc()    
    local g, FundRegionText, WaveParameterText, PhaseParameterText, PlotOptionsText, inputGenus, inputC, inputGraphicsChoice, inputMatrixParameters, inputPhaseParameters, inputPlotOptions, inputSolChoices,
    inputWaveParameters,RMPara,WavePara, c, PhasePara, U, U1, U2, V, Z, RM, SOLS, SOLSLIST, SOLSANUM, SOLSE, GOpt, num,PlotInfo,u;
    uses Maplets[Elements]:
    # Define the get input Maplet application. Define an abnormal
    # shutdown as one that returns FAIL.
    
    inputPlotOptions := Maplet( 'abnormalshutdown' = "FAIL",
        InputDialog['ID8'](
        PlotOptionsText,  
        'title' = "PlotOptions",
        'onapprove' = Shutdown( ['ID8'] ),
        'oncancel' = Shutdown( "FAIL" )
        )
    );
    
# Assign result to the value returned by the Maplet application.
    g := Maplets[Examples][GetInput]("Enter genus 1, 2 or 3:");
    if g = "FAIL" then
        error "Program Aborted";
    else
        try
            g := parse(g[1]);
            if g=1 then                
                FundRegionText := "Input a positive real number for the matrix";
                WaveParameterText := "Enter your choice of real numbers U, V";
                PhaseParameterText := "Enter your choice of real phase phi";
                                 
            elif g=2 then
                FundRegionText := "Input z11, z12, z22 so that 0<z11 <=z22, 0<=2z12<=z11";
                WaveParameterText := "Enter your choice of real numbers u1, u2 and v";
                PhaseParameterText := "Enter your choice of real phase phi1,phi2";
           
        
            elif g=3 then
                FundRegionText := "Input z11, z12, z13, z22, z23, z33 so that 0<z11<=z22<=z33, 0<2*z12<=z11, 0<2*z13<=z11, 2|z23|<=z22, 2(z12+z13-z23)<=z11+z22 ";
                WaveParameterText := "Enter your choice of real numbers u1, u2 and v";
                PhaseParameterText := "Enter your choice of real phase phi1, phi2, phi3";
           
            else
                error "Program Aborted: genus input is not 1, 2, or 3.";
            fi;
        catch:
            error "Program Aborted: genus input is not correct";
        end try; 
    fi;  
    print(cat("Genus:",convert(g,string)));
    RMPara:= Maplets[Examples][GetInput](FundRegionText);
    if RMPara = "FAIL" then
        error "Program Aborted";
    else
        RMPara := parse( RMPara );
        print(cat("Matrix Defined By:",convert([RMPara],string)));    
        WavePara := Maplets[Examples][GetInput](WaveParameterText);                      
        if WavePara = "FAIL" then
            error "Program Aborted";
        else            
            WavePara := parse( WavePara );
            print(cat("Wave Para defined By:",convert([WavePara],string)));  
            PhasePara:= Maplets[Examples][GetInput](PhaseParameterText);
            if PhasePara = "FAIL" then
                error "Program Aborted";
            else
                PhasePara :=  parse( PhasePara);
                print(cat("phase:",convert([PhasePara],string))); 
                c:= Maplets[Examples][GetInput]("Input a real number for c:");
                if c = "FAIL" then
                    error "Program Aborted";
                else        
                    c:=parse(c[1]);
                    print(cat("c:",convert(c,string))); 
                fi;
            fi;
        fi;
    fi;

    if g=1 then
        Z:= Matrix([RMPara]);
        B:=-Z/(2*Pi*I);
        SOLS:=waveConstants(B,WavePara):
    elif g=2 then
        Z:= Matrix(
            [
            [RMPara[1],RMPara[2]],
            [RMPara[2],RMPara[1]]
            ]
            )
            ;
        B:=-Z/(2*Pi*I);
        SOLS:=waveConstants(B,WavePara):
    else
        Z:= Matrix(
            [
            [RMPara[1],RMPara[2], RMPara[3]],
            [RMPara[2],RMPara[4], RMPara[5]],
            [RMPara[3],RMPara[5], RMPara[6]]
            ]
            )
            ;
        B:=-Z/(2*Pi*I);
        SOLS:=waveConstants(B,WavePara):
    fi;
    SOLSE:=SOLS[1][2];
    SOLSNUM:=SOLS[1][3];
    SOLSLIST:=SOLS[2..-1];
    SOLSANUM:=nops(SOLSLIST);
    for i from 1 to nops(SOLSLIST) do
        print(cat("solution ",convert(i,string)));
        print(SOLSLIST[i]);
        print("Error");
        print(SOLSE[i]);
    od; 
    GOpt:=Maplets[Examples][GetInput]("For plot type P, for Animation type A, otherwise type N:");
    print(cat("user chose",GOpt));
    if GOpt = "FAIL" then
            error "Program Aborted";
    fi;   
    if GOpt[1] = "P" then
        if SOLSANUM>1 then
            nums:= Maplets[Examples][GetInput]("From the list solutions printed, what solution(s) would you like to view? Enter you answer as a sequence of numbers."); 
            if nums = "FAIL" then
                error "Program Aborted";
            fi;
            nums:=parse(nums);        
            print(cat("Considering solutions", convert([nums],string)));  
        else 
            nums:=[1];
        fi; 
        PlotInfo := Maplets[Examples][GetInput]("Enter a,b,c,d,e, N where a<x<b, c<y<d, t=e, and N is the number of"); 
        if PlotInfo = "FAIL" then
            error "Program Aborted";
        fi; 
       PlotInfo := [parse(PlotInfo)];
       print(cat("PlotOptions", convert(PlotInfo,string)));  
            for int in nums do
                u := (x, y, t) -> KPgSol(x, y, t, B, SOLSLIST[int][1], SOLSLIST[int][2], SOLSLIST[int][3], [PhasePara], c);
                print(
                    plot3d(
                    Re(u(x, y, PlotInfo[5])),
                    x = PlotInfo[1] .. PlotInfo[2], y = PlotInfo[3] .. PlotInfo[4], colorscheme = ["Blue", "SkyBlue"], grid = [PlotInfo[6],PlotInfo[6]])
                    );
            od;
                    
    elif GOpt[1] = "A" then ;
        if SOLSANUM>1 then
            nums:=Maplets[Examples][GetInput]("From the list solutions printed, what solution would you like to animate."); 
            if nums = "FAIL" then
                error "Program Aborted";
            fi;  
            nums:=parse(nums);      
            print(cat("Considering solution", convert([nums],string)));  
        else
            nums:=1;
        fi; 
        PlotInfo := Maplets[Examples][GetInput]("Enter a,b,c,d,e,f, N, where a<x<b, c<y<d, e<t<f, and N is the number of "); 
        if PlotInfo = "FAIL" then
            error "Program Aborted";
        fi; 
       PlotInfo := [parse(PlotInfo)];             
       print(cat("PlotOptions", convert(PlotInfo,string)));      
        u := (x, y, t) -> KPgSol(x, y, t, B, SOLSLIST[nums][1], SOLSLIST[nums][2], SOLSLIST[nums][3], [PhasePara], c);
            print(
                plots[animate](
                plot3d,
                [Re(u(x, y, t)), x = PlotInfo[1] .. PlotInfo[2], y = PlotInfo[3] .. PlotInfo[4], colorscheme = ["Blue", "SkyBlue"], grid = [PlotInfo[7],PlotInfo[7]]],
                t= PlotInfo[5] .. PlotInfo[6] 
                )
                );            
    elif GOpt[1]="N" then
            
    return SOLS;
        
    else
        error "invalid choice for input";
      
    fi;
    return SOLS;
        



end proc;