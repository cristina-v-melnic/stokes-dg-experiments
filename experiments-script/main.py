import collections
import datetime
import os
import string
import numpy as np

###############################################################################
# variables to be adjusted, but which are the same for all tests

# string to preceed every produced input and output file
basename = "supg_param_opt"

# directory where all input files go
input_directory = "input_files"

# directory where all output files go (not the vtk files)
output_directory = "output_files"

# directory where all vtk files go
xdmf_directory = "xdmf_files"

# directory where the experiment results are stored
parent_dir = "../"

# file name of the executable
executable = "supg_param_opt"

# file name for pickle file for the parameter space
pickle_file_name = "parameter_space.pickle"

# directory where some postprocessing results can be written (e.g. plots)
results_dir = "results"

# file name for the file containing all relevant results as a pandas dataframe
h5_file_name = "data.h5"

# the date used for unique naming
date = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
###############################################################################

#_______________________________________________________________________________

# PART I: Preparation of the input ".dat" files needed for the simulations.

def create_dictionary(computation_dir):
    """
        Here we create a dictionary 'inputs' whose keys are strings which are
        ParMooN parameters. Their value is the value one wants to set in the
        resulting input file.
        
        :param computation_dir: (str) 
        	The directory where the results of this particular experiment will be stored.
        	Usually generated with the get_computation_dir() function, before calling 
        	create_dictionary() for the specific experiment.
    	
    	:return: (dict)
    		A dictionary of with all relevant parameters to add to the dat file in order
    		to run the simulation with parmoon. 
    """
    inputs = collections.OrderedDict()
    # Taken from nse_2d_dg_test.c++
    inputs["RE_NR"] = 10
    inputs["verbosity"] = 3
    inputs["example"] = 20
    inputs["problem_type"] = 3
    inputs["degree_polynomial"] = 0
    # only to check consistency
    inputs["is_RT_polynomial"] = "false"
    inputs["boundary_file"] = "Default_UnitSquare"
    inputs["geo_file"] = "TwoTriangles"
    inputs["symmetry_DG"] = 1
    inputs["face_sigma_DG"] = 10
    inputs["VELOCITY_SPACE"] = 1012
    inputs["refinement_n_initial_steps"] = 2
    inputs["outfile"] = parent_dir + \
        "{}/{}/default_parmoon_outfile.out".format(
            computation_dir, output_directory)
    inputs["output_directory"] = parent_dir + "/" + output_directory
    inputs["output_basename"] = "outfiles"

    inputs["separate_h5_files"] = "false"
    inputs["nonlinloop_maxit"] = 100
    inputs["nonlinloop_epsilon"] = 1e-10
    inputs["nonlinloop_slowfactor"] = 1.
    inputs["reynolds_number"] = 10.0

    inputs["output_write_xdmf"] = "true"
    inputs["output_xdmf_format"] = 1
    inputs["script_mode"] = "true"
    inputs["PRESSURE_SPACE"] = -4711
    inputs["NSTYPE"] = 3
    inputs["LAPLACETYPE"] = 0
    inputs["nse_nonlinear_form"] = "convective"
    inputs["space_discretization_type"] = "dg"
    inputs["output_compute_errors"] = "true"
    inputs["FLOW_PROBLEM_TYPE"] = 3
    inputs["USE_ISOPARAMETRIC"] = 0
    
    inputs["solver_verbosity"] = 2
    
    inputs["solver_type"] = "direct"
    inputs["direct_solver_type"] = "umfpack"

    #inputs["iterative_solver_type"] = "fgmres"

    inputs["residual_tolerance"] = 1.e-12
    #inputs["preconditioner"] = "vanka_nodal"
    inputs["INTERNAL_PROBLEM_LINEAR"] = 1

    #inputs["preconditioner"]= "multigrid"
    #inputs["multigrid_type"]= "mdml"
    #inputs["multigrid_n_levels"]= "1"
    #inputs["multigrid_smoother"] = "cell_vanka"
    #inputs["multigrid_smoother_coarse"] = "cell_vanka"
    inputs["max_n_iterations"] = 10000
    #inputs["multigrid_coarse_residual"] = 0.1
    #inputs["multigrid_vanka_damp_factor"] = 0.5
    #inputs["damping_factor"] = 0.5

    # added by me
    inputs["output_directory"] = parent_dir + \
        "{}/{}".format(computation_dir, xdmf_directory)

    return inputs


'''
def write_to_file(filename, dictionary):
    """ Write a given dictionary to a file.

    The file is newly created or overwritten. The dictionary can for
    example be a standard python dictionary or an OrderedDict from the
    collections module. The formatting is done such that ParMooN will be
    able to read this as an input file.
    """
    with open(filename, 'w') as inputfile:
        for param_name, param_value in dictionary.items():
            print("{}: {}".format(param_name, param_value), file=inputfile)
'''


def check_directory(dir):
    """ Check if directories exist and if not create them.
    
    :param dir: (str) 
    	The directory name to be created in the parent directory.
    """
    dir = os.path.join(parent_dir, dir)
    if not os.path.exists(dir):
        os.makedirs(dir)
    else:
        print("A directory {} already exists!".format(dir))
        exit(0)


def get_computation_dir(name = "None"):
    """
    Create a directory with the current time and date as part of the
    filename and return that name. 
    
    :param name: (str), optional 
        The name of the directory where the results of this particular
        experiment will be stored. By default it is given a name with 
        the date and the time of the run. 
       
    	
    :return:
    	The directory of the experiment containing other structured
    	directories needed for storing data, such as the input, output,
    	xdfm and results directories.
    	
    """

    # name of a directory where all output of this script goes
    if name != "None":
        computation_dir = name
    else:
        computation_dir = "computation_{}".format(date)

    check_directory(computation_dir)
    check_directory(computation_dir + "/" + input_directory)
    check_directory(computation_dir + "/" + output_directory)
    check_directory(computation_dir + "/" + xdmf_directory)
    check_directory(computation_dir + "/" + results_dir)
    return computation_dir


def find_keys_and_replace_values(file, keys, values, next_file):
    """
    This function copies the contents from "file" to "next_file" (can be different or the same,
    in the latter case, the content is overwritten), while changing the values of the parameters
    to correspond to the current dictionary from the method "create_dictionary()".

    :param file: (str)
    	The name of the file from which to copy the text and modify the values.
    	Usually one with the extension ".dat", found in the same directory as the main.py.
    	By default the "nse2d_dg_reference.dat" file must be used.
    :param keys: list of (str)
    	A list of keys corresponding to values that need to be changed in the desired 
    	simulation. If there is only one term, use the format [#key].
    :param values: list of data types appropriate to the parameters
    	A list of values that need to be changed in the new file, if there is only one
    	term, use the format [#value].
    :param next_file: 
    	The name of the file with changed values and the same text to be used to run
    	the desired parmoon simulation.
    
    :return:
    	A ".dat" file in the appropriate directory, ready for the desired parmoon simulation,
    	i.e., with changed selected parameters. 
    """

    with open(file, 'r') as f:
        text = f.readlines()

    f_next = open(next_file, "w")

    found = 0
    for line in text:
        for key in keys:
            if key + ":" in line:
                found += 1
                if (key + ": " + str(values[keys.index(key)]) + " ") in line:
                    f_next.write(line)
                else:
                    new_line = line.split(" ")
                    new_line[1] = str(values[keys.index(key)])
                    full_new_line = ""
                    for element in new_line:
                        full_new_line = full_new_line + element + " "
                    f_next.write(full_new_line+"\n")
        if found == 0:
            f_next.write(line)
        else:
            found = 0


def test_verbosity_values(initial_file):
    """ 
    Creates a new directory and 5 new files in input_files with
    different values of verbosity [1,5] and script_mode {true, false},
    all of them are copies of the initial_file with changed parameters.
    
    :param initial_file: (str)
    	The name of the reference ".dat" file containing the keys and values.
    
    :return:
    	A new directory and ".dat" files ready to be run for simulations
    	with different verbosities and adjustable script mode.
    """
    comp_dir = get_computation_dir()
    dictionary = create_dictionary(comp_dir)
    
    for i in range(1, 5):
        dictionary["verbosity"] = i
        find_keys_and_replace_values(initial_file, list(dictionary.keys()), list(dictionary.values()),
                                     comp_dir + "/" + input_directory + "/" + "nse2d_{}.dat".format(i))                             
        dictionary["script_mode"] = "false"

# END OF PART I: Input functions
#___________________________________________________________________________


#___________________________________________________________________________
# PART II: The simulation functions used for numerical experiments.

def get_refinement_experiment(directory_name = "None", example_nr = 20, velocity = 1011, RT = "false",
                                discret_type = "dg", nu_inv_factor = 10.0, face_sigma = 10, symmetry = 1):
    """
    Function to generate convergence histories, i.e, how the results of the 
    simulations evolve upon grid refinement.
    
    :param directory_name: (str)
    :param example_nr: (int)
    :param velocity: (int) 1___
    :param RT: {true, false}
    :param discret_type: (str)	
    :param nu_inv_factor: (float)
    :param face_sigma: (float)
    :param symmetry: (int) {-1, 0, 1}
    
    :return:
    	Parmoon simulation results for the given example ("example_nr"), using
    	the given velocity space ("velocity"), discretization type ("discret_type") 
    	and the sigma ("face_sigma") and "symmetry" for the case of the "dg" formulation.
    """
    comp_dir = get_computation_dir(directory_name)
    dictionary = create_dictionary(comp_dir)

    dictionary["VELOCITY_SPACE"] = velocity
    dictionary["is_RT_polynomial"] = RT
    dictionary["example"] = example_nr
    dictionary["space_discretization_type"] = discret_type
    dictionary["RE_NR"] = nu_inv_factor
    dictionary["reynolds_number"] = nu_inv_factor
    dictionary["face_sigma_DG"] = face_sigma
    dictionary["symmetry_DG"] = symmetry

    # if velocity == 101:
    #	dictionary["PRESSURE_SPACE"] = 1

    #refinement_n_initial_steps  [ 0, 20 ]
    for i in range(1,8):
        dictionary["refinement_n_initial_steps"] = i
        dictionary["output_basename"] = "ref_{}_sigma_{:.1e}_velocity_{}_example_{}".format(
            i, face_sigma, velocity, example_nr)
        dictionary["outfile"] = parent_dir + "{}/{}/".format(
            comp_dir, output_directory) + dictionary["output_basename"] + ".out"
        print("Experiment: {}, Space: {}, Refinement: {}, Re: {}".format(
            example_nr, velocity, i, nu_inv_factor))

        find_keys_and_replace_values("nse2d_dg_reference.dat", list(dictionary.keys()), list(dictionary.values()),
                                     parent_dir + "{}/{}/inputfile_{}.dat".format(comp_dir, input_directory, i))
        os.system("." + "/../" + "build/nse2d ../" + comp_dir +
                  "/" + input_directory + "/inputfile_{}.dat".format(i))


def loop_refinement_all_examples(experiment_name = "Refinement_Experiments"+"_{}".format(date)):
    """
    Function to generate a set of convergence histories with different
    examples, velocity spaces and discretization types in order to compare
    the optimality of different methods.
    
    :param experiment_name: (str), optional
    
    :return:
    	Parmoon simulation results in the following order experiment_directory,
    	a directory for each example, then a directory with each velocity space. 
    	
    """	
    check_directory(experiment_name)
    examples = collections.OrderedDict()

    examples = {}
    #examples[20] = "Polynomial.h"
    examples[2] = "SinCos.h"
    examples[5] = "polynomial_solution.h"

    for example_nr in examples.keys():
        example_name = examples[example_nr]
        check_directory(experiment_name + "/" + example_name)
        get_refinement_experiment(experiment_name + "/" + example_name +
                                  "/" + "BDM1", example_nr, velocity=1011, discret_type="dg")
        get_refinement_experiment(experiment_name + "/" + example_name +
                                  "/" + "RT0", example_nr, velocity=1000, discret_type="dg")
        get_refinement_experiment(experiment_name + "/" + example_name +
                                  "/" + "RT1", example_nr, velocity=1001, discret_type="dg")                          
        get_refinement_experiment(experiment_name + "/" + example_name +
                                  "/" + "MINI", example_nr, velocity=101, discret_type="galerkin")

        get_refinement_experiment(experiment_name + "/" + example_name +
                                  "/" + "BDM2", example_nr, velocity=1012, discret_type="dg")
        get_refinement_experiment(experiment_name + "/" + example_name +
                                  "/" + "RT2", example_nr, velocity=1002, discret_type="dg")
        get_refinement_experiment(experiment_name + "/" + example_name +
                                  "/" + "TH2", example_nr, velocity=2, discret_type="galerkin")

        get_refinement_experiment(experiment_name + "/" + example_name +
                                  "/" + "BDM3", example_nr, velocity=1013, discret_type="dg")
        get_refinement_experiment(experiment_name + "/" + example_name +
                                  "/" + "RT3", example_nr, velocity=1003, discret_type="dg")
        get_refinement_experiment(experiment_name + "/" + example_name +
                                  "/" + "TH3", example_nr, velocity=3, discret_type="galerkin")

        get_refinement_experiment(experiment_name + "/" + example_name + "/" + "SV2", example_nr, velocity = 12, discret_type = "galerkin")
        get_refinement_experiment(experiment_name + "/" + example_name + "/" + "SV3", example_nr, velocity = 13, discret_type = "galerkin")
        get_refinement_experiment(experiment_name + "/" + example_name + "/" + "SV4", example_nr, velocity = 14, discret_type = "galerkin")


def get_Re_nr_scaling(experiment_name = "Reynolds_nr_Scaling"+"_{}".format(date)):
    """
    Function to study pressure robustness by generating sets of convergence
    histories with different Reynolds number and observing whether or not they 
    differ.
    
    :param experiment_name: (str), optional
    
    :return:
    	Parmoon simulation results in the following order experiment_directory,
    	a directory for each example, then a directory for every combination of
    	a particular Re and each velocity space. 
    	
    """	
    check_directory(experiment_name)
    examples = collections.OrderedDict()

    examples = {}
    #examples[20] = "Polynomial.h"
    examples[2] = "SinCos.h"
    examples[5] = "polynomial_solution.h"

    nu_inv_factors = [1.0, 100.0, 10000.0, 1000000.0]

    for example_nr in examples.keys():
        
        example_name = examples[example_nr]
        check_directory(experiment_name + "/" + example_name)
        
        for nu_inv in nu_inv_factors:
            get_refinement_experiment(experiment_name + "/" + example_name + "/" + "BDM1_Re={}".format(
                nu_inv), example_nr, velocity=1011, discret_type="dg", nu_inv_factor=nu_inv)
            get_refinement_experiment(experiment_name + "/" + example_name + "/" + "MINI_Re={}".format(
                nu_inv), example_nr, velocity=101, discret_type="galerkin", nu_inv_factor=nu_inv)
            get_refinement_experiment(experiment_name + "/" + example_name + "/" + "RT0_Re={}".format(
                nu_inv), example_nr, velocity=1000, discret_type="dg", nu_inv_factor=nu_inv)
            get_refinement_experiment(experiment_name + "/" + example_name + "/" + "RT1_Re={}".format(
                nu_inv), example_nr, velocity=1001, discret_type="dg", nu_inv_factor=nu_inv)    

            get_refinement_experiment(experiment_name + "/" + example_name + "/" + "BDM2_Re={}".format(
                nu_inv), example_nr, velocity=1012, discret_type="dg", nu_inv_factor=nu_inv)
            get_refinement_experiment(experiment_name + "/" + example_name + "/" + "TH2_Re={}".format(
                nu_inv), example_nr, velocity=2, discret_type="galerkin", nu_inv_factor=nu_inv)
            get_refinement_experiment(experiment_name + "/" + example_name + "/" + "RT2_Re={}".format(
                nu_inv), example_nr, velocity=1002, discret_type="dg", nu_inv_factor=nu_inv)  
                
            get_refinement_experiment(experiment_name + "/" + example_name + "/" + "BDM3_Re={}".format(
                nu_inv), example_nr, velocity=1013, discret_type="dg", nu_inv_factor=nu_inv)
            get_refinement_experiment(experiment_name + "/" + example_name + "/" + "TH3_Re={}".format(
                nu_inv), example_nr, velocity=3, discret_type="galerkin", nu_inv_factor=nu_inv)
            get_refinement_experiment(experiment_name + "/" + example_name + "/" + "RT3_Re={}".format(
                nu_inv), example_nr, velocity=1003, discret_type="dg", nu_inv_factor=nu_inv)       
                
            get_refinement_experiment(experiment_name + "/" + example_name + "/" + "SV2_Re={}".format(
                nu_inv), example_nr, velocity = 12, discret_type = "galerkin", nu_inv_factor=nu_inv)
            get_refinement_experiment(experiment_name + "/" + example_name + "/" + "SV3_Re={}".format(
                nu_inv), example_nr, velocity = 13, discret_type = "galerkin", nu_inv_factor=nu_inv)


def check_optimal_sigma(RE_NR = 10):
    """
    Function to determine the optimal choice of the sigma parameter for dg 
    by comparing the magnitudes of errors from simulations across ranges of 
    values for sigma.
    	    
    :param RE_NR: (float), optional
    
    :return:
    	Parmoon simulation results in the following order: 3 experiment directories,
    	for each range of sigma values, a directory for each example, then a directory
    	for every combination of a particular sigma value and each velocity space.
    """		
    
    examples = collections.OrderedDict()

    examples = {}
    #examples[20] = "Polynomial.h"
    examples[2] = "SinCos.h"
    examples[5] = "polynomial_solution.h"

    # for the first search
    ps = np.arange(-7, 8)
    sigma = np.power(np.ones(len(ps))*RE_NR, ps)
    
    loop_sigma_search(sigma, examples, experiment_name = "Optimal_sigma_ref1"+"_{}".format(date))


    # for the 2nd search
    multiples = np.arange(1, 11)
    sigma = np.multiply(np.ones(len(multiples))*RE_NR, multiples)
    
    loop_sigma_search(sigma, examples, experiment_name = "Optimal_sigma_ref2"+"_{}".format(date))


    # for the 3rd search
    multiples = np.arange(5, 36, 5)
    sigma = np.multiply(np.ones(len(multiples)), multiples)
    
    loop_sigma_search(sigma, examples, experiment_name = "Optimal_sigma_ref3"+"_{}".format(date))

    

def loop_sigma_search(sigma, examples, experiment_name):
    """
    Function to generate convergence histories for every value
    of sigma used in the function check_optimal_sigma().
    
    :param RE_NR: list of (float)
    :param examples: list of (str)
    :param experiment_name: (str)

    :return:
    	Simulation results for all sigma values in the directory
    	structure given by the get_refinement_experiment() function.
    """
    check_directory(experiment_name)
	
    for example_nr in examples.keys():
        example_name = examples[example_nr]
    check_directory(experiment_name + "/" + example_name)
    for s in range(len(sigma)):
            # BDM1 with all types of symmetry SIPG, NIPG and IPG
        get_refinement_experiment(experiment_name + "/" + example_name + "/" + "BDM1_sigma={}_symm={}".format(
                sigma[s], 1), example_nr, velocity = 1011, discret_type = "dg", face_sigma = sigma[s], symmetry = 1)
        get_refinement_experiment(experiment_name + "/" + example_name + "/" + "BDM1_sigma={}_symm={}".format(
                sigma[s], -1), example_nr, velocity = 1011, discret_type = "dg", face_sigma = sigma[s], symmetry = -1)
        get_refinement_experiment(experiment_name + "/" + example_name + "/" + "BDM1_sigma={}_symm={}".format(
                sigma[s], 0), example_nr, velocity = 1011, discret_type = "dg", face_sigma = sigma[s], symmetry = 0)
            
            # RT1 with all types of symmetry SIPG, NIPG and IPG
        get_refinement_experiment(experiment_name + "/" + example_name + "/" + "RT1_sigma={}_symm={}".format(
                sigma[s], 1), example_nr, velocity = 1001, discret_type = "dg", face_sigma = sigma[s], symmetry = 1)
        get_refinement_experiment(experiment_name + "/" + example_name + "/" + "RT1_sigma={}_symm={}".format(
                sigma[s], -1), example_nr, velocity = 1001, discret_type = "dg", face_sigma = sigma[s], symmetry = -1)
        get_refinement_experiment(experiment_name + "/" + example_name + "/" + "RT1_sigma={}_symm={}".format(
                sigma[s], 0), example_nr, velocity = 1001, discret_type = "dg", face_sigma = sigma[s], symmetry = 0)	
	

# END OF PART II: Experiment functions
#_____________________________________________________________________________________________________________	

if __name__ == '__main__':
    
    # To obtain convergence history data for specified examples, uncomment below. 	
    loop_refinement_all_examples()
    
    # To obtain convergence history data for different Re values to check for 
    # pressure robustness uncomment below.
    #get_Re_nr_scaling()   	
    
    # To obtain convergence history data for different values of sigma to find
    # the appropriate parameter uncomment below.
    #check_optimal_sigma()
    
    
