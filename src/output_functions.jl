using PrettyTables,Plots,Statistics,Parameters

"""
        describe(gl::ModelGL,gl_2::ModelGL)

Print solution status as well as values and description of parameters and aggregate values:
Comparison of initial and terminal steady-state.
"""
function describe(gl::ModelGL,gl_2::ModelGL)
    print_status(gl,gl_2 ) 
    print_params(gl,gl_2 ) 
end

"""
        describe(gl::ModelGL)   

Print solution status as well as values and description of parameters and aggregate values.
"""
function describe(gl::ModelGL)
    print_status(gl ) 
    print_params(gl ) 
end

"""
        print_params(gl::ModelGL)

Print values and description of parameters and aggregate values.
"""
function print_params(gl::ModelGL) 

        params = ["β","ν", "ϕ","γ","η","fin","sep"]
        desc   = ["discount","UI benefits", "borr. constraint", 
                        "risk aversion", "eta", "job-finding rate", 
                        "job-separation rate", ] 
        vals   = [ getfield(gl,Symbol(i) )   for i in params ]
 
        data = hcat((params,vals,desc)... )
 

        tgt_nam  = ["NE","νY", "B_4Y","D1_4Y","D2_4Y"]
        tgt_dsc  = ["avg. hours of employed","UI benefit / GDP",
                        "liquid wealth / GDP", "HH debt-to-GDP (initial ss)",
                        "HH debt-to-GDP (terminal ss)" ]
        tgt_vals = [ getfield(gl,Symbol(i) )   for i in tgt_nam ]

        if gl.steadystate_type == "initial"
                actual_nam  = ["NE_actual","νY_actual", "B_4Y_actual","D_4Y_actual" ] 
                actual_vals = [[ getfield(gl,Symbol(i) )   for i in actual_nam ]...; NaN ]
        elseif  gl.steadystate_type == "terminal"
                actual_nam  = ["NE_actual","νY_actual", "B_4Y_actual", "D_4Y_actual"] 
                actual_vals = [[ getfield(gl,Symbol(i) )   for i in actual_nam[1:3] ]...; NaN; getfield(gl,Symbol(actual_nam[end])) ]
        end

        params = ["β","ν", "ϕ","γ","η","fin","sep","ρ","σϵ"]
        desc   = ["Discount factor","Unemployment benefits", "Borrowing constraint", 
                        "Risk aversion", "Leisure curvature", "Job-finding rate", 
                        "Job-separation rate","Persistance income","SD income shocks"  ]
        vals   = [ getfield(gl,Symbol(i) )   for i in params ]
        data = hcat((params,vals ,fill("",length(params)),desc)... )
 

        tgt_nam  = ["NE","νY", "B_4Y","D1_4Y","D2_4Y" ]
        tgt_dsc  = ["Average hours of employed","Unemployment benefits per GDP",
                "Liquid wealth per GDP", "Household debt-to-GDP (initial ss)",
                "Household debt-to-GDP (terminal ss)"  ]
        tgt_vals =  [ getfield(gl,Symbol(i) )   for i in tgt_nam  ] 
        actual_nam      = ["NE_actual","νY_actual", "B_4Y_actual","D_4Y_actual" ] 

        if (gl.steadystate_type == "initial")   
                actual_vals  = [[ getfield(gl,Symbol(i) )   for i in actual_nam ]...; NaN;   ]
        elseif (gl.steadystate_type == "terminal") 
                actual_vals  =  [[ getfield(gl ,Symbol(i) )   for i in actual_nam[1:3]]...; NaN; getfield(gl ,Symbol(actual_nam[end]))  ]
        end

         
        data2 = hcat((tgt_nam,actual_vals ,tgt_vals,tgt_dsc)... )
 
        filler1 = ["Parameters" ""  "" "" ]
        filler2 = [" " ""   "" ""; "Aggregates" ""   "" ""]
        pretty_table(vcat((filler1,data,filler2,data2)...), tf = tf_simple,  
        formatters = ft_printf("%5.3f") ,alignment=[:l,:r,:r, :l],
        header = (["Field", "Value",  "Target","Description" ],
                  [" ", gl.steadystate_type , " "," " ] ),
        highlighters = ( hl_cell( [(1,1);(length(params)+3,1)], crayon"bold") ),
        border_crayon = crayon"bold yellow",  header_crayon = crayon"bold green",  subheader_crayon = crayon"green",
        crop = :none)

end

"""
        print_params(gl::ModelGL,gl_2::ModelGL)

Print values and description of parameters and aggregate values: 
Comparison of initial and terminal steady-state.
"""
function print_params(gl::ModelGL,gl_2::ModelGL) 

        params = ["β","ν", "ϕ","γ","η","fin","sep","ρ","σϵ"]
        desc   = ["Discount factor","Unemployment benefits", "Borrowing constraint", 
                "Risk aversion", "Leisure curvature", "Job-finding rate", 
                "Job-separation rate","Persistance income","SD income shocks" ]
        vals_1   = [ getfield(gl,Symbol(i) )   for i in params ]
        vals_2   = [ getfield(gl_2,Symbol(i) )   for i in params ]

        data = hcat((params,vals_1,vals_2,fill("",length(params)),desc)... )
 

        tgt_nam  = ["NE","νY", "B_4Y","D1_4Y","D2_4Y" ]
        tgt_dsc  = ["Average hours of employed","Unemployment benefits per GDP",
                        "Liquid wealth per GDP", "Household debt-to-GDP (initial ss)",
                        "Household debt-to-GDP (terminal ss)"  ]
        tgt_vals =  [ getfield(gl,Symbol(i) )   for i in tgt_nam  ] 

        actual_nam      = ["NE_actual","νY_actual", "B_4Y_actual","D_4Y_actual" ] 
        actual_vals_1 = [[ getfield(gl,Symbol(i) )   for i in actual_nam ]...; NaN;   ]
        actual_vals_2 = [[ getfield(gl_2,Symbol(i) )   for i in actual_nam[1:3]]...; NaN; getfield(gl_2,Symbol(actual_nam[end]))  ]
 
        data2 = hcat((tgt_nam,actual_vals_1,actual_vals_2,tgt_vals,tgt_dsc)... )
 
        filler1 = ["Parameters" "" ""  "" ""]
        filler2 = [" " "" ""  "" ""; "Aggregates" "" ""  "" ""]
        pretty_table(vcat((filler1,data,filler2,data2)...), tf = tf_simple,  
        formatters = ft_printf("%5.3f"),alignment=[:l,:r,:r,:r,:l],
        header = (["Field", "Value 1)", "Value 2)","Target","Description" ],
                  [" ", gl.steadystate_type, gl_2.steadystate_type, "(1)", "" ] ),
        highlighters = ( hl_cell( [(1,1);(length(params)+3,1)], crayon"bold") ),

        border_crayon = crayon"yellow",  header_crayon = crayon"bold green",  subheader_crayon = crayon"green",
        crop = :none)

end
 
"""
        print_status(gl::ModelGL) 

Print solution status.
"""
function print_status(gl::ModelGL) 


        data = ["Model" ""   "";
                "Steady-state solved    ?" "."^10 "$(gl.steadystate_solved)"; 
                "Calibrated             ?" "."^10 "$(gl.calibrated)";
                "Steady-state type      ?" "."^10 "$(gl.steadystate_type)";
                "" ""   "" ]
 
        pretty_table(data, tf = tf_borderless,
                        noheader = true,
                        cell_alignment = Dict( (1,1) => :l, (2,1) => :l, (3,1) => :l, (4,1) => :l ),
                        formatters = ft_printf("%10.1f", 2),
                        highlighters = (hl_cell( [(1,1);(6,1)], crayon"bold"),
                                        hl_col(2, crayon"dark_gray")),
                        body_hlines = [1,6],
                        body_hlines_format = Tuple('─' for _ = 1:4) )
                   

end


"""
        print_status(gl::ModelGL,gl_2::ModelGL) 

Print solution status of both models.
"""
function print_status(gl::ModelGL,gl_2::ModelGL) 


        data = ["1) Model" ""   "";
                "Steady-state solved    ?" "."^10 "$(gl.steadystate_solved)"; 
                "Calibrated             ?" "."^10 "$(gl.calibrated)";
                "Steady-state type      ?" "."^10 "$(gl.steadystate_type)";
                "" ""   "" ]
        data2 = ["2) Model" ""   "";
                "Steady-state solved    ?" "."^10 "$(gl_2.steadystate_solved)"; 
                "Calibrated             ?" "."^10 "$(gl_2.calibrated)";
                "Steady-state type      ?" "."^10 "$(gl_2.steadystate_type)";
                "" ""   "" ]
        data = hcat(data,data2)
        pretty_table(data, tf = tf_borderless,
                        noheader = true,
                        cell_alignment = Dict( (1,1) => :l, (2,1) => :l, (3,1) => :l, (4,1) => :l ,
                                               (1,4) => :l, (2,4) => :l, (3,4) => :l, (4,4) => :l   ),
                        formatters = ft_printf("%10.1f", 2),
                        highlighters = (hl_cell( [(1,1);(1,4)], crayon"bold"),
                                        hl_col(2, crayon"dark_gray")),
                        body_hlines = [1,6],
                        body_hlines_format = Tuple('─' for _ = 1:4) )
                       

end




"""
        (F1,p1,p2) = plots_figure_1(gl::ModelGL)

Plots figure I graphs.        
"""
function plots_figure_1(gl::ModelGL)
        plotly(linewidth = 2.)
  
        @unpack b_grid,pr = gl                  # objects that are identical across both steady states
 
        # Figure I (from paper)

        Y1 = gl.Y_actual
        p1 = plot( b_grid[b_grid .>= -gl.ϕ]/(4*Y1), gl.c_pol[2,b_grid .>= -gl.ϕ],label = "θ = 2")
        plot!(b_grid[b_grid .>= -gl.ϕ]/(4*Y1), gl.c_pol[8,b_grid .>= -gl.ϕ],
                linestyle = :dash,label = "θ = 8")
        title!("consumption")
        xaxis!([-gl.ϕ,12.5])
 

        p2 = plot( b_grid[b_grid .>= -gl.ϕ]/(4*Y1), gl.n_pol[2,gl.b_grid .>= -gl.ϕ],label =  "θ = 2")
        plot!(b_grid[b_grid .>= -gl.ϕ]/(4*Y1), gl.n_pol[8,gl.b_grid .>= -gl.ϕ],
                linestyle = :dash,label = "θ = 8")
        title!("labor supply")
        xaxis!([-gl.ϕ,12.5])
 
        F1 = plot(p1, p2, layout = (1, 2), legend = false) 
        display(F1)
        return (F1,p1,p2)
end


"""
        (F4,p1,p2) = plots_figure_4(gl::ModelGL,gl_tss::ModelGL)

Plots figure IV graphs.        
"""
function plots_figure_4(gl::ModelGL,gl_tss::ModelGL)
        plotly(linewidth = 2.)
  
        @unpack b_grid,pr = gl                  # objects that are identical across both steady states
 
        # Figure IV (from paper)
        Y1 = gl.Y_actual
        b_acc1 =  ((pr'*gl.b_pol[:,b_grid     .>= -gl.ϕ])'     .- b_grid[b_grid .>= -gl.ϕ])/(4*Y1)
        b_acc2 =  ((pr'*gl_tss.b_pol[:,b_grid .>= -gl.ϕ])' .- b_grid[b_grid .>= -gl.ϕ])/(4*Y1)
        p1 =plot( b_grid[b_grid .>= -gl.ϕ]/(4*Y1), b_acc1)
        plot!(b_grid[b_grid .>= -gl.ϕ]/(4*Y1), b_acc2,
                linestyle = :dash)
        title!("bond accumulation policy")
        xaxis!([-2,14],-2:2:14)
        yaxis!([-0.4,0.6],-0.4:0.2:0.6)
 

        bond_distribution1 = dropdims( pr'    *gl.JD[:,b_grid .>= -gl.ϕ]     ,dims=1)
        bond_distribution2 = dropdims( pr'*gl_tss.JD[:,b_grid .>= -gl.ϕ] ,dims=1)
        p2 =plot( b_grid[b_grid .>= -gl.ϕ]/(4*Y1),  bond_distribution1,label = "initial")
        plot!(b_grid[b_grid .>= -gl.ϕ]/(4*Y1),  bond_distribution2,label = "terminal",
                linestyle = :dash)
        title!("bond distribution")
        xaxis!([-2,14],-2:2:14)
        yaxis!([0,0.004])
        xaxis!([-gl.ϕ,12.5])
 
        F4 = plot(p1, p2, layout = (2, 1), legend = false) 
        display(F4)

        return (F4,p1,p2)
end


"""
        (F3,p1,p2,p3,p4)  = plots_figure_4(gl::ModelGL,gl_tss::ModelGL)

Plots figure III graphs.        
"""
function plots_figure_3(gl::ModelGL,Tgl::TransGL)
        plotly(linewidth = 2.)
        Y1 = gl.Y_actual
 
        # Figure III (from paper)
        Tp = 24    # number of periods plotted

        p1 = plot(0:Tp, Tgl.ϕ_t[1:Tp+1]./(4*Tgl.Y_t[1:Tp+1]),ylims = (0.5,1), yticks = 0.5:0.1:10, legend = false, title = "borrowing limit")         # borrowing limit
        p2 = plot(0:Tp, Tgl.D_4Y_t[1:Tp+1], ylims = (0.08,0.2), yticks = 0:0.02:0.2, legend = false, title = "household debt-to-GDP ratio")             # debt2gdp ratio
        p3 = plot(0:Tp, [gl.r;Tgl.r_t[1:Tp]].*400, ylims = (-2,2.5), yticks = -2:0.5:2, legend = false, title = "interest rate")       # annualized interest rate
        p4 = plot(0:Tp, [0, 100*(Tgl.Y_t[1:Tp+1]./Y1.-1)], ylims = (-1.2,0), yticks = -1:0.2:0, legend = false, title = "output")   # output deviation from steady state
        F3 = plot(p1, p2, p3, p4, layout = (2, 2), legend = false)
 
        display(F3)

        return (F3,p1,p2,p3,p4) 
end

















 
