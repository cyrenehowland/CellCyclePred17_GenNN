/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2022, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <ctime>
#include <cmath>
#include <omp.h>
#include <fstream>

#include "./core/PhysiCell.h"
#include "./modules/PhysiCell_standard_modules.h" 

// put custom code modules here! 

#include "./custom_modules/custom.h" 
	
using namespace BioFVM;
using namespace PhysiCell;

#include <cmath>  // For exp(), pow()


void add_random_substrate_to_voxel()
{
    int num_voxels = microenvironment.number_of_voxels();
    int random_voxel_index = rand() % num_voxels;
    
    int idx_food = microenvironment.find_density_index("food");
    
    // Gaussian distribution parameters
    double mean_x = microenvironment.mesh.voxels[random_voxel_index].center[0];
    double mean_y = microenvironment.mesh.voxels[random_voxel_index].center[1];
//    double mean_z = microenvironment.mesh.voxels[random_voxel_index].center[2];
    double sigma = 100.0 + (rand() % 10); // Standard deviation in microns
    double amplitude = 200.0 + (rand() % 10); // Amplitude of the Gaussian distribution
    
    // Apply Gaussian distribution to a single randomly chosen voxel
    for( int n = 0; n < num_voxels; n++ )
    {
        std::vector<double> voxel_center = microenvironment.mesh.voxels[n].center;
        double x = voxel_center[0];
        double y = voxel_center[1];
//        double z = voxel_center[2];
        
        // Calculate the Gaussian function value based on the distance from the random voxel
//        double distance_squared = pow(x - mean_x, 2) + pow(y - mean_y, 2) + pow(z - mean_z, 2);
        double distance_squared = pow(x - mean_x, 2) + pow(y - mean_y, 2);
        double value = amplitude * exp(-distance_squared / (2 * pow(sigma, 2)));
        
        // Add Gaussian value to the existing substrate concentration
        microenvironment(n)[idx_food] += value;  // Correct way to access and modify substrate value
    }
}

//
//void add_linear_gradient_substrate()
//{
//    // Find the index for the new nutrient source (e.g., "food")
//    int idx_food = microenvironment.find_density_index("food");
//
//    // Get the size of the microenvironment mesh
//    int num_voxels = microenvironment.number_of_voxels();
//
//    // Define the maximum concentration at the corner (0,0,0)
//    double max_concentration = 200.0;  // Adjust this value as needed
//
//    // Define the distance scale for the gradient (larger values make the gradient more gradual)
//    double distance_scale = 500.0;  // Adjust to control how quickly the gradient decreases
//
//    // Loop through all the voxels and set a gradient from the origin
//    for (int n = 0; n < num_voxels; n++)
//    {
//        // Get the center coordinates of the voxel
//        std::vector<double> voxel_center = microenvironment.mesh.voxels[n].center;
//        double x = voxel_center[0];
//        double y = voxel_center[1];
//        double z = voxel_center[2];
//
//        // Calculate distance from the origin (0,0,0)
//        double distance_from_origin = sqrt(x * x + y * y + z * z);
//
//        // Calculate the concentration based on the linear gradient
//        double concentration = max_concentration * (1.0 - (distance_from_origin / distance_scale));
//
//        // Ensure the concentration is non-negative
//        if (concentration < 0.0)
//        {
//            concentration = 0.0;
//        }
//
//        // Update the microenvironment with the calculated concentration using a Dirichlet node
//        microenvironment.update_dirichlet_node(n, idx_food, concentration);
//    }
//}

//void toggle_sunlight_gradient(double current_time, double sunlight_on_duration, double sunlight_off_duration)
//{
//    // Find the index for the "sunlight" substrate
//    int idx_sunlight = microenvironment.find_density_index("food");
//
//    // Calculate the current cycle duration (on + off period)
//    double cycle_duration = sunlight_on_duration + sunlight_off_duration;
//
//    // Determine the phase of the cycle (on or off)
//    double time_in_cycle = fmod(current_time, cycle_duration);
//
//    // If time_in_cycle is within the 'on' duration, apply the gradient
//    if (time_in_cycle < sunlight_on_duration)
//    {
//        // Define the maximum concentration at the corner (0,0,0)
//        double max_concentration = 4.5;  // Adjust this value as needed
//
//        // Define the distance scale for the gradient (larger values make the gradient more gradual)
//        double distance_scale = 1000.0;  // Adjust to control how quickly the gradient decreases
//
//        // Apply the gradient as before
//        int num_voxels = microenvironment.number_of_voxels();
//        for (int n = 0; n < num_voxels; n++)
//        {
//            std::vector<double> voxel_center = microenvironment.mesh.voxels[n].center;
//            double x = voxel_center[0];
//            double y = voxel_center[1];
//            double z = voxel_center[2];
//
//            // Calculate distance from the origin (0,0,0)
//            double distance_from_origin = sqrt(x * x + y * y + z * z);
//
//            // Calculate the concentration based on the linear gradient
//            double concentration = max_concentration * (1.0 - (distance_from_origin / distance_scale));
//
//            // Ensure the concentration is non-negative
//            if (concentration < 0.0)
//            {
//                concentration = 0.0;
//            }
//
//            // Update the microenvironment with the calculated concentration using a Dirichlet node
//            microenvironment.update_dirichlet_node(n, idx_sunlight, concentration);
//        }
//    }
//    else
//    {
//
//        // If time_in_cycle is in the 'off' phase, set the sunlight concentration to zero
//        int num_voxels = microenvironment.number_of_voxels();
//        for (int n = 0; n < num_voxels; n++)
//        {
//            microenvironment.update_dirichlet_node(n, idx_sunlight, 0.0);
//        }
//    }
//}

void toggle_sunlight_gradient(double current_time, double sunlight_on_duration, double sunlight_off_duration)
{
    // Find the index for the "sunlight" substrate
    int idx_sunlight = microenvironment.find_density_index("food");

    // Calculate the current cycle duration (on + off period)
    double cycle_duration = sunlight_on_duration + sunlight_off_duration;

    // Determine the phase of the cycle (on or off)
    double time_in_cycle = fmod(current_time, cycle_duration);

    // If time_in_cycle is within the 'on' duration, apply the gradient
    if (time_in_cycle < sunlight_on_duration)
    {
        // Define the maximum concentration at the -x wall
        double max_concentration = 3.0;  // Adjust this value as needed

        // Define the distance scale for the gradient (larger values make the gradient more gradual)
        double distance_scale = 2000.0;  // Adjust to control how quickly the gradient decreases

        // Apply the gradient
        int num_voxels = microenvironment.number_of_voxels();
        for (int n = 0; n < num_voxels; n++)
        {
            std::vector<double> voxel_center = microenvironment.mesh.voxels[n].center;
            double x = voxel_center[0];  // Use the x coordinate for gradient calculation

            // Calculate the distance from the -x axis (e.g., from a wall at the negative x-boundary)
            double distance_from_x_wall = fabs(x + distance_scale);  // Adjust if needed based on your coordinate system

            // Calculate the concentration based on the linear gradient
            double concentration = max_concentration * (1.0 - (distance_from_x_wall / distance_scale));

            // Ensure the concentration is non-negative
            if (concentration < 0.0)
            {
                concentration = 0.0;
            }

            // Update the microenvironment with the calculated concentration using a Dirichlet node
            microenvironment.update_dirichlet_node(n, idx_sunlight, concentration);
        }
    }
    else
    {
        // If time_in_cycle is in the 'off' phase, set the sunlight concentration to zero
        int num_voxels = microenvironment.number_of_voxels();
        for (int n = 0; n < num_voxels; n++)
        {
            microenvironment.update_dirichlet_node(n, idx_sunlight, 0.0);
        }
    }
}


int main( int argc, char* argv[] )
{
	// load and parse settings file(s)
	
	bool XML_status = false; 
	char copy_command [1024]; 
	if( argc > 1 )
	{
		XML_status = load_PhysiCell_config_file( argv[1] ); 
		sprintf( copy_command , "cp %s %s" , argv[1] , PhysiCell_settings.folder.c_str() ); 
	}
	else
	{
		XML_status = load_PhysiCell_config_file( "./config/PhysiCell_settings.xml" );
		sprintf( copy_command , "cp ./config/PhysiCell_settings.xml %s" , PhysiCell_settings.folder.c_str() ); 
	}
	if( !XML_status )
	{ exit(-1); }
	
	// copy config file to output directry 
	system( copy_command ); 
	
	// OpenMP setup
	omp_set_num_threads(PhysiCell_settings.omp_num_threads);
	
	// time setup 
	std::string time_units = "min"; 

	/* Microenvironment setup */ 
	
	setup_microenvironment(); // modify this in the custom code 
	
	/* PhysiCell setup */ 
 	
	// set mechanics voxel size, and match the data structure to BioFVM
	double mechanics_voxel_size = 30; 
	Cell_Container* cell_container = create_cell_container_for_microenvironment( microenvironment, mechanics_voxel_size );
	
	/* Users typically start modifying here. START USERMODS */ 
	
	create_cell_types();
	
	setup_tissue();

	/* Users typically stop modifying here. END USERMODS */ 
	
	// set MultiCellDS save options 

	set_save_biofvm_mesh_as_matlab( true ); 
	set_save_biofvm_data_as_matlab( true ); 
	set_save_biofvm_cell_data( true ); 
	set_save_biofvm_cell_data_as_custom_matlab( true );
	
	// save a simulation snapshot 
	
	char filename[1024];
	sprintf( filename , "%s/initial" , PhysiCell_settings.folder.c_str() ); 
	save_PhysiCell_to_MultiCellDS_v2( filename , microenvironment , PhysiCell_globals.current_time ); 
	
	// save a quick SVG cross section through z = 0, after setting its 
	// length bar to 200 microns 

	PhysiCell_SVG_options.length_bar = 200; 

	// for simplicity, set a pathology coloring function 
	
	std::vector<std::string> (*cell_coloring_function)(Cell*) = my_coloring_function; 
	
	sprintf( filename , "%s/initial.svg" , PhysiCell_settings.folder.c_str() ); 
	SVG_plot( filename , microenvironment, 0.0 , PhysiCell_globals.current_time, cell_coloring_function );
	
	sprintf( filename , "%s/legend.svg" , PhysiCell_settings.folder.c_str() ); 
	create_plot_legend( filename , cell_coloring_function ); 
	
	display_citations(); 
	
	// set the performance timers 

	BioFVM::RUNTIME_TIC();
	BioFVM::TIC();
	
	std::ofstream report_file;
	if( PhysiCell_settings.enable_legacy_saves == true )
	{	
		sprintf( filename , "%s/simulation_report.txt" , PhysiCell_settings.folder.c_str() ); 
		
		report_file.open(filename); 	// create the data log file 
		report_file<<"simulated time\tnum cells\tnum division\tnum death\twall time"<<std::endl;
	}
	
	// main loop 
	
	try 
	{		
		while( PhysiCell_globals.current_time < PhysiCell_settings.max_time + 0.1*diffusion_dt )
		{
			// save data if it's time. 
			if( fabs( PhysiCell_globals.current_time - PhysiCell_globals.next_full_save_time ) < 5.0 * diffusion_dt )
			{
				display_simulation_status( std::cout ); 
				if( PhysiCell_settings.enable_legacy_saves == true )
				{	
					log_output( PhysiCell_globals.current_time , PhysiCell_globals.full_output_index, microenvironment, report_file);
				}
				
				if( PhysiCell_settings.enable_full_saves == true )
				{	
					sprintf( filename , "%s/output%08u" , PhysiCell_settings.folder.c_str(),  PhysiCell_globals.full_output_index ); 
					
					save_PhysiCell_to_MultiCellDS_v2( filename , microenvironment , PhysiCell_globals.current_time ); 
				}
				
				PhysiCell_globals.full_output_index++; 
				PhysiCell_globals.next_full_save_time += PhysiCell_settings.full_save_interval;
			}
			
			// save SVG plot if it's time
			if( fabs( PhysiCell_globals.current_time - PhysiCell_globals.next_SVG_save_time  ) < 5.0 * diffusion_dt )
			{
				if( PhysiCell_settings.enable_SVG_saves == true )
				{	
					sprintf( filename , "%s/snapshot%08u.svg" , PhysiCell_settings.folder.c_str() , PhysiCell_globals.SVG_output_index ); 
					SVG_plot( filename , microenvironment, 0.0 , PhysiCell_globals.current_time, cell_coloring_function );
					
					PhysiCell_globals.SVG_output_index++; 
					PhysiCell_globals.next_SVG_save_time  += PhysiCell_settings.SVG_save_interval;
				}
			}

			// update the microenvironment
			microenvironment.simulate_diffusion_decay( diffusion_dt );
            
            toggle_sunlight_gradient(PhysiCell_globals.current_time, 600, 600);

			
//            // Add food randomly into arena periodically
//            if( fabs(fmod(PhysiCell_globals.current_time, 200.0)) < 1e-2 )
//            {
//                // Periodically add substrate to voxels
//                add_random_substrate_to_voxel();
////                add_linear_gradient_substrate();
////                toggle_sunlight_gradient(PhysiCell_globals.current_time, 100, 200);
//                std::cout << std::endl << "Substrate Added " << std::endl;
//            }
            
			// run PhysiCell
			((Cell_Container *)microenvironment.agent_container)->update_all_cells( PhysiCell_globals.current_time );
			
			/*
			  Custom add-ons could potentially go here. 
			*/
			
			PhysiCell_globals.current_time += diffusion_dt;
		}
		
		if( PhysiCell_settings.enable_legacy_saves == true )
		{			
			log_output(PhysiCell_globals.current_time, PhysiCell_globals.full_output_index, microenvironment, report_file);
			report_file.close();
		}
	}
	catch( const std::exception& e )
	{ // reference to the base of a polymorphic object
		std::cout << e.what(); // information from length_error printed
	}
	
	// save a final simulation snapshot 
	
	sprintf( filename , "%s/final" , PhysiCell_settings.folder.c_str() ); 
	save_PhysiCell_to_MultiCellDS_v2( filename , microenvironment , PhysiCell_globals.current_time ); 
	
	sprintf( filename , "%s/final.svg" , PhysiCell_settings.folder.c_str() ); 
	SVG_plot( filename , microenvironment, 0.0 , PhysiCell_globals.current_time, cell_coloring_function );
	
	// timer 
	
	std::cout << std::endl << "Total simulation runtime: " << std::endl; 
	BioFVM::display_stopwatch_value( std::cout , BioFVM::runtime_stopwatch_value() ); 

	return 0; 
}
