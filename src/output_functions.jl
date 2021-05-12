using PrettyTables


# ===========  WORK IN PROGRESS  ===========

function describe(gl::ModelGL)

    print_status(gl ) 
    print_params(gl ) 

end


function print_params(gl::ModelGL) 

    params = ["β","ν", "ϕ","γ","η","fin","sep"]
    desc   = ["discount","UI benefits", "borr. constraint", 
                "risk aversion", "eta", "job-finding rate", "job-separation rate"]
    vals   = [ getfield(gl,Symbol(i) )   for i in params ]
    data = hcat(hcat(hcat(params,vals),vals),desc)
    
    pretty_table(data, tf = tf_simple, title = "Parameters", formatters = ft_printf("%5.3f"),
                header = (["Parameter", "Value (initial)", "Value (terminal)","Description"] ),
                border_crayon = crayon"bold yellow",  header_crayon = crayon"bold green")

 
end

function print_status(gl::ModelGL) 
    
    data = ["Status" ""   "";
            "Policy solved?" "."^10 "yes/no"; 
            "Calibrated?" "."^10 "not solved/default/calibrated";
            "Steady-state?" "."^10 "not solved/initial/terminal";
            "" ""   ""
            "Other" ""  "" ]
    
    pretty_table(data, tf = tf_borderless,
                noheader = true,
                cell_alignment = Dict( (1,1) => :l, (6,1) => :l ),
                formatters = ft_printf("%10.1f", 2),
                highlighters = (hl_cell( [(1,1);(6,1)], crayon"bold"),
                                hl_col(2, crayon"dark_gray")),
                body_hlines = [1,6],
                body_hlines_format = Tuple('─' for _ = 1:4) )
 
end