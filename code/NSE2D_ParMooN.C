#include <Database.h>
#include <FEDatabase2D.h>
#include <ParameterDatabase.h>
#include <Chrono.h>


#include <iostream>
#include "NavierStokes_DG.h"
#include "NavierStokes.h"

//for the date
#include <stdio.h>
#include <time.h>

// for output to_string() funct
#include <string>


int main(int, char* argv[])
{
  Chrono timer; // start a stopwatch which measures time spent in program parts
  
  //  declaration of database, you need this in every program
  TDatabase::create(argv[1]);
  TFEDatabase2D FEDatabase;
  ParameterDatabase parmoon_db = ParameterDatabase::parmoon_default_database();
  parmoon_db.read(argv[1]);
  	
  
  //open OUTFILE, this is where all output is written to (additionally to console)
  Output::set_outfile(parmoon_db["outfile"], parmoon_db["script_mode"]);
  Output::setVerbosity(parmoon_db["verbosity"]);
  
  // default construct a domain object
  TDomain domain(parmoon_db);
  
  // possibly change parameters in the database, if they are not meaningful now
  check_parameters_consistency_NSE(parmoon_db);
  
  // write all Parameters to the OUTFILE (not to console) for later reference
  parmoon_db.write(Output::get_outfile());
  TDatabase::WriteParamDB(argv[0]);

  // refine grid
  domain.refine_and_get_hierarchy_of_collections(parmoon_db);
  
  // write grid into an Postscript file
  if(parmoon_db["output_write_ps"])
    domain.PS("Domain.ps", It_Finest, 0);
    
  // create an object of the Navier-Stokes class with dg or galerkin spatial discretization
  if(parmoon_db["space_discretization_type"].is("dg"))
  {
  	NavierStokes_DG<2> nse2d(domain, parmoon_db);
	nse2d.assemble();
	nse2d.solve();
	
        timer.restart_and_print("solving procedure");
	// if parameter_db["output_basename"]. doesn't work
	// std::string name = parmoon_db["output_directory"].get<std::string>()+ "/" + parmoon_db["output_basename"].get<std::string>()+"_min_max.out"; 
        
        
        std::vector<std::shared_ptr<FEMatrix>> blocks = nse2d.get_matrix().get_blocks_uniquely();
	for (int i=0; i<4; i++)
	{
		Output::redirect("../../../software/matlab_program/bin/files/search_1/matrix_part_" + std::to_string(i) + "_" + parmoon_db["output_basename"].get<std::string>());
		blocks.at(i).get()->PrintFull();
		Output::resetOutfile();
	}
	nse2d.output();
        	
  }
  else if(parmoon_db["space_discretization_type"].is("galerkin"))
  {
       NavierStokes<2> nse2d(domain, parmoon_db);
       nse2d.assemble_linear_terms();
       nse2d.solve();
       timer.restart_and_print("solving procedure");
       
       std::vector<std::shared_ptr<FEMatrix>> blocks = nse2d.get_matrix().get_blocks_uniquely();
	for (int i=0; i<4; i++)
	{
		Output::redirect("../../../software/matlab_program/bin/files/search_1/matrix_part_" + std::to_string(i) + "_" + parmoon_db["output_basename"].get<std::string>());
		blocks.at(i).get()->PrintFull();
		Output::resetOutfile();
	}
       nse2d.output();	
  }
  else
  {
  	cout << "This space_discretization_type is not supported yet.";
  	return 0;
  }
  
  timer.print_total_time("complete NSE2D_ParMooN program");
  TDatabase::destroy();
  Output::close_file();
  return 0;
}
