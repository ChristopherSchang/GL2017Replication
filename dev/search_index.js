var documenterSearchIndex = {"docs":
[{"location":"functions/#Functions","page":"Functions","title":"Functions","text":"","category":"section"},{"location":"functions/","page":"Functions","title":"Functions","text":"Modules = [GL2017Replication]","category":"page"},{"location":"functions/#GL2017Replication.ModelGL","page":"Functions","title":"GL2017Replication.ModelGL","text":"Structure holding all parameters and solutions to the steady-state model\n\n\n\n\n\n","category":"type"},{"location":"functions/#GL2017Replication.TransGL","page":"Functions","title":"GL2017Replication.TransGL","text":"Structure holding transition objects\n\n\n\n\n\n","category":"type"},{"location":"functions/#GL2017Replication.EGM!-Tuple{ModelGL}","page":"Functions","title":"GL2017Replication.EGM!","text":"EGM!(gl::ModelGL)\n\nCalculates optimal policy functions with the EGM algorithm.\n\n\n\n\n\n","category":"method"},{"location":"functions/#GL2017Replication.EGM_trans!-Tuple{ModelGL, TransGL, Int64}","page":"Functions","title":"GL2017Replication.EGM_trans!","text":"EGM_trans!(gl::ModelGL,Tgl::TransGL,t::Int)\n\nCalculates optimal policy functions with one step of the EGM algorithm.\n\n\n\n\n\n","category":"method"},{"location":"functions/#GL2017Replication.aggregate!-Tuple{ModelGL}","page":"Functions","title":"GL2017Replication.aggregate!","text":"aggregate!(gl::ModelGL;terminal::Bool = false)\n\nComputes aggregate values and their distance from target. Requires that the joint distribution has already been solved (compute_distribution!) If ´terminal = true´ only the terminal debt target will be considered.\n\n\n\n\n\n","category":"method"},{"location":"functions/#GL2017Replication.calibrate!-Tuple{ModelGL}","page":"Functions","title":"GL2017Replication.calibrate!","text":"calibrate!(gl::ModelGL)\n\nCalibrate the model to the target values. Does not require any prior function calls.\n\n\n\n\n\n","category":"method"},{"location":"functions/#GL2017Replication.calibrate_terminal-Tuple{ModelGL}","page":"Functions","title":"GL2017Replication.calibrate_terminal","text":"gl_tss = calibrate_terminal(gl_initial::ModelGL)\n\nCalibrate the model to the terminal target value (only borrowing constraint) Requires a calibrated initital steady-state object. Returns a new terminal steady-state object.\n\n\n\n\n\n","category":"method"},{"location":"functions/#GL2017Replication.compute_distribution!-Tuple{ModelGL}","page":"Functions","title":"GL2017Replication.compute_distribution!","text":"compute_distribution!(gl::ModelGL)\n\nComputes the joint distribution of productivity and bond holdings. Requires that optimal policies have already been solved (EGM!)\n\n\n\n\n\n","category":"method"},{"location":"functions/#GL2017Replication.compute_steady_state!-Tuple{ModelGL}","page":"Functions","title":"GL2017Replication.compute_steady_state!","text":"compute_steady_state!(gl::ModelGL)\n\nComputes steady-state policy functions with given set of parameters.\n\n\n\n\n\n","category":"method"},{"location":"functions/#GL2017Replication.describe-Tuple{ModelGL, ModelGL}","page":"Functions","title":"GL2017Replication.describe","text":"    describe(gl::ModelGL,gl_2::ModelGL)\n\nPrint solution status as well as values and description of parameters and aggregate values: Comparison of initial and terminal steady-state.\n\n\n\n\n\n","category":"method"},{"location":"functions/#GL2017Replication.describe-Tuple{ModelGL}","page":"Functions","title":"GL2017Replication.describe","text":"    describe(gl::ModelGL)\n\nPrint solution status as well as values and description of parameters and aggregate values.\n\n\n\n\n\n","category":"method"},{"location":"functions/#GL2017Replication.domath-Tuple{Number}","page":"Functions","title":"GL2017Replication.domath","text":"domath(x::Number)\n\nReturn x + 5.\n\n\n\n\n\n","category":"method"},{"location":"functions/#GL2017Replication.find_cl-NTuple{9, Any}","page":"Functions","title":"GL2017Replication.find_cl","text":"F = find_cl(c, j, b1, b2, r, θ, z, fac, gameta)\n\nHelper function to find consumption at the constraint. Returns the residual the budget constraint.\n\nc - consumption\nj - state index\nb1 - bond holdings today\nb2 - bond holdings tomorrow\nr - interest rate\nθ - state value array\nfac - parameter\ngameta - parameter\n\nReturn\n\nF - residual budget constraint\n\n\n\n\n\n","category":"method"},{"location":"functions/#GL2017Replication.hello-Tuple{String}","page":"Functions","title":"GL2017Replication.hello","text":"hello(who::String)\n\nReturn \"Hello, who\".\n\n\n\n\n\n","category":"method"},{"location":"functions/#GL2017Replication.initilize!-Tuple{ModelGL}","page":"Functions","title":"GL2017Replication.initilize!","text":"initilize!(gl::ModelGL)\n\nAdds borrowing constraint as gridpoint\n\n\n\n\n\n","category":"method"},{"location":"functions/#GL2017Replication.load_parameters_iss!-Tuple{Any}","page":"Functions","title":"GL2017Replication.load_parameters_iss!","text":"load_parameters_iss!(gl))\n\nPre-load calibrated parameters to speed up calibration.\n\n\n\n\n\n","category":"method"},{"location":"functions/#GL2017Replication.load_parameters_tss!-Tuple{Any}","page":"Functions","title":"GL2017Replication.load_parameters_tss!","text":"load_parameters_tss!(gl))\n\nPre-load calibrated parameters to speed up calibration (terminal).\n\n\n\n\n\n","category":"method"},{"location":"functions/#GL2017Replication.print_params-Tuple{ModelGL, ModelGL}","page":"Functions","title":"GL2017Replication.print_params","text":"    print_params(gl::ModelGL,gl_2::ModelGL)\n\nPrint values and description of parameters and aggregate values:  Comparison of initial and terminal steady-state.\n\n\n\n\n\n","category":"method"},{"location":"functions/#GL2017Replication.print_params-Tuple{ModelGL}","page":"Functions","title":"GL2017Replication.print_params","text":"    print_params(gl::ModelGL)\n\nPrint values and description of parameters and aggregate values.\n\n\n\n\n\n","category":"method"},{"location":"functions/#GL2017Replication.print_status-Tuple{ModelGL, ModelGL}","page":"Functions","title":"GL2017Replication.print_status","text":"    print_status(gl::ModelGL,gl_2::ModelGL)\n\nPrint solution status of both models.\n\n\n\n\n\n","category":"method"},{"location":"functions/#GL2017Replication.print_status-Tuple{ModelGL}","page":"Functions","title":"GL2017Replication.print_status","text":"    print_status(gl::ModelGL)\n\nPrint solution status.\n\n\n\n\n\n","category":"method"},{"location":"functions/#GL2017Replication.transition!-Tuple{ModelGL, ModelGL, TransGL}","page":"Functions","title":"GL2017Replication.transition!","text":"transition!(gl::ModelGL,gl_tss::ModelGL,Tgl::TransGL)\n\nComputes transition path.\n\n\n\n\n\n","category":"method"},{"location":"functions/","page":"Functions","title":"Functions","text":"end","category":"page"},{"location":"transition/#Transition","page":"Transition","title":"Transition","text":"","category":"section"},{"location":"transition/","page":"Transition","title":"Transition","text":"This part is of particular interest. When we tried to replicate Figure 3 of the original paper without loading the transition objects previously obtained by the authors (i.e., we set rerun neq 1), but keeping the default values of speed = 05 and decay = 03 (two parameters governing the updating rule of ), we couldn't make the code converge. It seems to us that the authors used two different parameters that are a priori impossible to deduce. Importantly, for one of the few combinations of speed and decay that would produce output (i.e., 0.1 and 0.01) in the original MATLAB file, we obtain essentially the same evolution for all variables included in Figure 3 (trivially, the borrowing constraint, but also household debt, the interest rate, and aggregate output). This somehow reassures us about the validity of our approach, although we cannot obviously rule out completely the possibility that we missed some other features of the code explaining the puzzle.  Below, we show our reproduction of Figure 3 in the original paper. ","category":"page"},{"location":"transition/","page":"Transition","title":"Transition","text":"end","category":"page"},{"location":"steadystate/#Steady-States","page":"Steady State","title":"Steady States","text":"","category":"section"},{"location":"steadystate/#Solving-the-Model","page":"Steady State","title":"Solving the Model","text":"","category":"section"},{"location":"steadystate/","page":"Steady State","title":"Steady State","text":"An easy way to solve the model in steady-state with the given set of default parameters is shown below. Note that all the steps have to occur sequentially. First, instantiate the ModelGL structure that holds all parameter values and solutions of the model.","category":"page"},{"location":"steadystate/","page":"Steady State","title":"Steady State","text":"gl = ModelGL() ","category":"page"},{"location":"steadystate/","page":"Steady State","title":"Steady State","text":"Calling ","category":"page"},{"location":"steadystate/","page":"Steady State","title":"Steady State","text":"compute_steady_state!(gl)","category":"page"},{"location":"steadystate/","page":"Steady State","title":"Steady State","text":"will solve for the steady-state policy function, joint distribution over productivity and asset states as well as aggregate variables.","category":"page"},{"location":"steadystate/#Calibrating-the-Steady-State","page":"Steady State","title":"Calibrating the Steady-State","text":"","category":"section"},{"location":"steadystate/#Initial-Steady-State","page":"Steady State","title":"Initial Steady-State","text":"","category":"section"},{"location":"steadystate/","page":"Steady State","title":"Steady State","text":"To calibrate the model in steady-state to a given set of target values do the following: First, instantiate the ModelGL structure as above.","category":"page"},{"location":"steadystate/","page":"Steady State","title":"Steady State","text":"gl = ModelGL() ","category":"page"},{"location":"steadystate/","page":"Steady State","title":"Steady State","text":"Calling the calibrate function will solve the model many times and calibrate the parameters to the target values:","category":"page"},{"location":"steadystate/","page":"Steady State","title":"Steady State","text":"calibrate!(gl)","category":"page"},{"location":"steadystate/","page":"Steady State","title":"Steady State","text":"The solution and the calibrated parameters will be stored in the ModelGL structure.","category":"page"},{"location":"steadystate/#Terminal-Steady-State","page":"Steady State","title":"Terminal Steady-State","text":"","category":"section"},{"location":"steadystate/","page":"Steady State","title":"Steady State","text":"To calibrate the model to the terminal steady-state we first need a calibrated model for the initial steady-state (step above):","category":"page"},{"location":"steadystate/","page":"Steady State","title":"Steady State","text":"gl = ModelGL()  \ncalibrate!(gl)","category":"page"},{"location":"steadystate/","page":"Steady State","title":"Steady State","text":"Given the solved initial steady-state we can calibrate the terminal steady-state to the new debt-target","category":"page"},{"location":"steadystate/","page":"Steady State","title":"Steady State","text":"gl_tss = calibrate_terminal(gl)","category":"page"},{"location":"steadystate/","page":"Steady State","title":"Steady State","text":"where gl_tss denotes the terminal steady-state object.","category":"page"},{"location":"steadystate/#Results","page":"Steady State","title":"Results","text":"","category":"section"},{"location":"steadystate/","page":"Steady State","title":"Steady State","text":"At any point we can query basic statistics (parameter values, aggregates and targets) with the describe function","category":"page"},{"location":"steadystate/","page":"Steady State","title":"Steady State","text":"describe(gl)","category":"page"},{"location":"steadystate/","page":"Steady State","title":"Steady State","text":"or","category":"page"},{"location":"steadystate/","page":"Steady State","title":"Steady State","text":"describe(gl,gl_tss)","category":"page"},{"location":"steadystate/","page":"Steady State","title":"Steady State","text":"for a side-by-side comparison of the initial and terminal steady-state.","category":"page"},{"location":"steadystate/","page":"Steady State","title":"Steady State","text":"end","category":"page"},{"location":"#Replication-of-Guerrieri,-Lorenzoni-(QJE,2017)","page":"Home","title":"Replication of Guerrieri, Lorenzoni (QJE,2017)","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"This replication study was part of our evaluation for the course Numerical Methods at SciencesPo Paris in Spring 2021","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"#About-The-Project","page":"Home","title":"About The Project","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"This replication study replicates parts of the paper ´Credit crises, precautionary savings, and the liquidity trap´ (Guerrieri, Lorenzoni (2017)). The code notation mostly follows the original version of the authors in MATLAB.","category":"page"},{"location":"#Structure","page":"Home","title":"Structure","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The Package is build around two structures: ","category":"page"},{"location":"","page":"Home","title":"Home","text":"The ModelGL structure that holds all parameters and (if solved/calibrated) the policy functions and aggregate variables.","category":"page"},{"location":"","page":"Home","title":"Home","text":"gl = ModelGL() ","category":"page"},{"location":"","page":"Home","title":"Home","text":"At any point you query the status of the model with","category":"page"},{"location":"","page":"Home","title":"Home","text":"describe(gl)","category":"page"},{"location":"","page":"Home","title":"Home","text":"The TransGL structure that holds the objects for the transition from an initital to a terminal steady-state","category":"page"},{"location":"","page":"Home","title":"Home","text":"gl_trans = TransGL() ","category":"page"},{"location":"#References","page":"Home","title":"References","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Guerrieri, V. & Lorenzoni, G. Credit crises, precautionary savings, and the liquidity trap  The Quarterly Journal of Economics, Oxford University Press, 2017, 132, 1427-1467","category":"page"},{"location":"","page":"Home","title":"Home","text":"end","category":"page"}]
}
