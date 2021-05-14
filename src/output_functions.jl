using PrettyTables
 
function describe(gl::ModelGL,gl_2::ModelGL)
    print_status(gl,gl_2 ) 
    print_params(gl,gl_2 ) 
end

function describe(gl::ModelGL)
    print_status(gl ) 
    print_params(gl ) 
end

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
 

function print_status(gl::ModelGL ) 


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




 
