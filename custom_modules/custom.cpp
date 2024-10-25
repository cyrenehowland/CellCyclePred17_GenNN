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
# Copyright (c) 2015-2021, Paul Macklin and the PhysiCell Project             #
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
#include "./custom_intracellular.h"
#include "./custom.h"
#include <cstdlib>
#include <cmath>
#include <iostream>






void create_cell_types( void )
{
    // set the random seed
    SeedRandom( parameters.ints("random_seed") );
    
    /*
       Put any modifications to default cell definition here if you
       want to have "inherited" by other cell types.
       
       This is a good place to set default functions.
    */
    
    initialize_default_cell_definition();
    cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment );
    
    cell_defaults.functions.volume_update_function = NULL;
    cell_defaults.functions.update_velocity = standard_update_cell_velocity;

    cell_defaults.functions.update_migration_bias = NULL;
    cell_defaults.functions.update_phenotype = NULL; // update_cell_and_death_parameters_O2_based;
    cell_defaults.functions.custom_cell_rule = NULL;
    cell_defaults.functions.contact_function = NULL;
    
    cell_defaults.functions.add_cell_basement_membrane_interactions = NULL;
    cell_defaults.functions.calculate_distance_to_membrane = NULL;
    
    /*
       This parses the cell definitions in the XML config file.
    */
    
    initialize_cell_definitions_from_pugixml();

    /*
       This builds the map of cell definitions and summarizes the setup.
    */
        
    build_cell_definitions_maps();

    /*
       This intializes cell signal and response dictionaries
    */

    setup_signal_behavior_dictionaries();

    /*
       Cell rule definitions
    */

    setup_cell_rules();

    /*
       Put any modifications to individual cell definitions here.
       
       This is a good place to set custom functions.
    */
    
    
//    cell_defaults.functions.update_phenotype = phenotype_function;
//    cell_defaults.functions.custom_cell_rule = NULL;
//    cell_defaults.functions.contact_function = contact_function;
    
     // Make sure we're ready for 2D
    cell_defaults.functions.set_orientation = up_orientation;
    cell_defaults.phenotype.geometry.polarity = 1.0;
    cell_defaults.phenotype.motility.restrict_to_2D = true;
    
    
    

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // ~~~~~~~~~~~~~ (Prey) Cell Function Assignment ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Cell_Definition* pPreyDef = find_cell_definition( "Cell" );
    pPreyDef->functions.update_phenotype = prey_phenotype_function;
//    pPreyDef->phenotype.mechanics.attachment_elastic_constant = parameters.doubles("attachment_elastic_constant");
    pPreyDef->functions.contact_function = standard_elastic_contact_function;
//    pPreyDef -> functions.volume_update_function = prey_growth_and_metabolism;
    pPreyDef -> functions.cell_division_function = custom_division_function;


    // Cell Cycle:
    Cycle_Model& prey_cycle_model = pPreyDef -> phenotype.cycle.model(); // Access cycle model for modification

    // Don't let division, entry functions or exit functions occur automatically after any phase
    for (int i = 0; i < prey_cycle_model.phases.size(); ++i)
    {
        prey_cycle_model.phases[i].division_at_phase_exit = false; // reset division flag for all phases
        prey_cycle_model.phases[i].entry_function = NULL;

    }
  

    
//    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//    // ~~~~~~~~~~~~~ Predator Cell Function Assignment ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Cell_Definition* pPredDef = find_cell_definition( "Predator" );

    Cycle_Model& pred_cycle_model = pPredDef -> phenotype.cycle.model(); // Access cycle model for modification
    
    // Don't let division occur automatically after any phase
    for (int i = 0; i < pred_cycle_model.phases.size(); ++i)
    {
        pred_cycle_model.phases[i].division_at_phase_exit = false; // reset division flag for all phases
    }

    pPredDef->functions.update_phenotype = pred_phenotype_function;
    pPredDef -> functions.volume_update_function = standard_volume_update_function;

    
    /*
       This builds the map of cell definitions and summarizes the setup.
    */
    display_cell_definitions( std::cout );
    
    return;
}


void prey_initialize_properties(Cell* cell){
    cell -> custom_data["initial_volume"] = cell -> phenotype.volume.total;
    cell -> custom_data["energy"] = 25;
    cell -> phenotype.motility.is_motile = true;
//    cell -> custom_data["pred_stress"] = 0;
    cell -> custom_data["danger"] = 0;
//    cell -> custom_data["attachment_loss_time"]= -1.0;
    
//    std::cout << "Cell " << cell->ID << " , INITIALIZED AND attachment loss time  = " << cell->custom_data["attachment_loss_time"] << std::endl;
    
    
}

//void pred_initialize_properties(Cell* cell){
//    cell -> custom_data["initial_volume"] = cell -> phenotype.volume.total;
//    cell -> custom_data["energy"] = 10;
//    cell -> phenotype.motility.is_motile = true;
//    cell -> custom_data["edible"] = -1.0;
//
//}



void initialize_cell_network(Cell* cell, int input_size, int hidden_size, int output_size) {
    // unique random seed: ensures unique initialization for each cell
    srand(cell->ID);
    // create and asign neural network
    cell->phenotype.intracellular = new NeuralNetworkIntracellular(input_size, hidden_size, output_size);
    // Use dynamic_cast to check the cast
    auto* nn = dynamic_cast<NeuralNetworkIntracellular*>(cell->phenotype.intracellular);
    
    // Check if the cast was successful
    if (nn == nullptr) {
        std::cerr << "Error: Failed to cast phenotype.intracellular to NeuralNetworkIntracellular." << std::endl;
        return;
    }
    nn->initialize_random(cell->ID);
    nn -> print_weights(); // Print weights after initialization
    
    // Initialize custom_data for inputs and outputs
    for (size_t i = 0; i < input_size; ++i) {
        std::string input_var_name = "input_" + std::to_string(i);
        cell->custom_data.add_variable(input_var_name, "input " + std::to_string(i), 0.0);
        std::cout << "Added custom data variable: " << input_var_name << " to cell " << cell->ID << std::endl;

    }

    for (size_t i = 0; i < output_size; ++i) {
        std::string output_var_name = "output_" + std::to_string(i);
        cell->custom_data.add_variable(output_var_name, "output " + std::to_string(i), 0.0);
        std::cout << "Added custom data variable: " << output_var_name << " to cell " << cell->ID << std::endl;

    }
    
}

std::vector<double> get_cell_inputs(Cell* cell) {
    
    // Smooth out/prepare inputs:
    
    // Input indicator if cell is currently attached
    bool attached = false;
    if (cell->state.number_of_attached_cells() > 0){attached = true;}
    
    // Parameters for Hill function for danger
//    double half_max_danger = 10.0;  // Half-maximal concentration
//    double hill_pow_danger = 2.0;  // Hill power
    // Apply Hill function to calculate pred_stress from danger level
//    double danger = cell->custom_data["danger"];
//    std::cout << "Cell ID (" << cell->ID << ") has danger = " << danger << std::endl;
//    cell->custom_data["pred_stress"] = Hill_response_function(danger, half_max_danger, hill_pow_danger);
//    std::cout << "Cell ID (" << cell->ID << ") has pred stress = " << cell->custom_data["pred_stress"] << std::endl;

    
    // Parameters for Hill function for energy
    double half_max_energy = 100.0;  // Half-maximal concentration
    double hill_pow_energy = 2.0;  // Hill power
    // Apply Hill function to energy
    double hill_energy = Hill_response_function(cell->custom_data["energy"], half_max_energy, hill_pow_energy);
//    std::cout << "Cell ID (" << cell->ID << ") has hill energy = " << hill_energy << std::endl;
    
    // Input for nutrient availible nearby
    double nutrient = cell->nearest_density_vector()[0];
//    std::cout << "Cell ID (" << cell->ID << ") has nutrient = " << nutrient << std::endl;
    
    // Have input boolean for if there is any nutrient detected
    bool day = false;
    if (nutrient > 0){day = true;}
    else {day = false;}
    
    // Make volume into a multiple of initial_volume
    double initial_volume = cell -> custom_data["initial_volume"];
    double volume_multiple = cell->phenotype.volume.total / initial_volume;
    
    // Make age relative to age killed in death condtion
    double age = (cell -> custom_data["time_since_last_division"])/5000;
    
    // Reserve space for 7 inputs to avoid reallocations
    std::vector<double> inputs;
    inputs.reserve(6);  // You know the size ahead of time
    
    // Now add the inputs
    inputs.push_back(nutrient);
    inputs.push_back(day);
    inputs.push_back(volume_multiple);
    inputs.push_back(hill_energy);
    inputs.push_back(age);
//    inputs.push_back(cell->custom_data["pred_stress"]);
    inputs.push_back(attached);
    
    
    

//    inputs.push_back(cell -> custom_data["time_spent_attached"]);
    
    return inputs;
}



void setup_microenvironment( void )
{
	// set domain parameters 
	
	// put any custom code to set non-homogeneous initial conditions or 
	// extra Dirichlet nodes here. 
	
	// initialize BioFVM 
	
	initialize_microenvironment(); 	
	
	return; 
}

void setup_tissue( void )
{
	double Xmin = microenvironment.mesh.bounding_box[0]; 
	double Ymin = microenvironment.mesh.bounding_box[1]; 
	double Zmin = microenvironment.mesh.bounding_box[2]; 

	double Xmax = microenvironment.mesh.bounding_box[3]; 
	double Ymax = microenvironment.mesh.bounding_box[4]; 
	double Zmax = microenvironment.mesh.bounding_box[5]; 
	
	if( default_microenvironment_options.simulate_2D == true )
	{
		Zmin = 0.0; 
		Zmax = 0.0; 
	}
	
	double Xrange = Xmax - Xmin; 
	double Yrange = Ymax - Ymin; 
	double Zrange = Zmax - Zmin; 
	
	// create some of each type of cell 
	
	Cell* pC;
	
	for( int k=0; k < cell_definitions_by_index.size() ; k++ )
	{
		Cell_Definition* pCD = cell_definitions_by_index[k]; 
		std::cout << "Placing cells of type " << pCD->name << " ... " << std::endl; 
		for( int n = 0 ; n < parameters.ints("number_of_cells") ; n++ )
		{
			std::vector<double> position = {0,0,0}; 
			position[0] = Xmin + UniformRandom()*Xrange; 
			position[1] = Ymin + UniformRandom()*Yrange; 
			position[2] = Zmin + UniformRandom()*Zrange; 
			
			pC = create_cell( *pCD ); 
			pC->assign_position( position );
		}
	}
	std::cout << std::endl; 
    
    // load cells from your CSV file (if enabled)
    load_cells_from_pugixml();
    
    
    // Initialize cell Properties:
    for( auto& cell : *all_cells )
    {
        // Check if the current cell is a default cell
        if (cell->type == find_cell_definition("Cell")->type){
            
            // Initilaize cell properties:
            prey_initialize_properties(cell);
            initialize_cell_network((cell), 6, 15, 2); // Example sizes: 3 inputs, 5 hidden neurons, 2 outputs
        }
    
        
//        if (cell->type == find_cell_definition("Predator")->type){
//            
//            // Initilaize cell properties:
//            pred_initialize_properties(cell);
//        }
    }
	
	return; 
}



double time_since_last_division(Cell* pCell)
{
    double last_division_time = pCell->custom_data["last_division_time"];
    double current_time = PhysiCell_globals.current_time;

    return current_time - last_division_time;
}



std::vector<std::string> my_coloring_function( Cell* pCell )
{ return paint_by_number_cell_coloring(pCell); }


void evaluate_death_conditions(Cell* pCell, Phenotype& phenotype, double dt){
    
    // Death based on starvation:
    if( pCell->custom_data["energy"] < 1.0 )
    {
        pCell->lyse_cell();
        return;
    }
    
    
    // Death based on size constraint (if too small or too big, die)
    double initial_volume = pCell->custom_data["initial_volume"];
    double current_volume = pCell->phenotype.volume.total;
    if (current_volume < initial_volume || current_volume > 30 * initial_volume)
    {
        pCell->lyse_cell();
        return;
    }
    
    
    // Death based on time since last division:
    if (pCell -> custom_data["time_since_last_division"] > 5000){
        //        std::cout << "Cell ID (" << pCell->ID << ") Will die now. Its been too long since last division" << std::endl;
        pCell->lyse_cell();
        return;
        
    }
    
    // Define the probability of lysis for unattached cells
    double lysis_probability = 0.009; // Adjust this value to tune how likely unattached cells are to lyse
    
    // Check if the cell is attached
    if (pCell->state.number_of_attached_cells() == 0) // Unattached
    {
        // Generate a random number between 0 and 1
        double random_number = UniformRandom();
        
        // If the random number is less than the lysis probability, lyse the cell
        if (random_number < lysis_probability)
        {
            pCell->lyse_cell();
        }
    }
}

       // Additional logic for attached cells, if needed...
   
    
    
//    // Get the index for apoptosis
//    int apoptosis_index = pCell->phenotype.death.find_death_model_index( "Apoptosis" );
//
//    // Check if the cell is attached
//    if (pCell->state.number_of_attached_cells() == 0) // Unattached
//    {
//        // Increase death rate for unattached cells
//        phenotype.death.rates[apoptosis_index] = 0.0005;
//    }
//    else // Attached
//    {
//        // Normal death rate for attached cells
//        phenotype.death.rates[apoptosis_index] = 0;
//    }





void prey_phenotype_function( Cell* pCell, Phenotype& phenotype, double dt )
{
    
    // Calculate time since last division
    pCell->custom_data["time_since_last_division"] = time_since_last_division(pCell);
    
    if (pCell->state.number_of_attached_cells() > 0){
        // If you are attached to others, label 'inedible'
        pCell -> custom_data["edible"] = false;
    }
    else {
        pCell -> custom_data["edible"] = true;
    }
    

    // Neural Network processing
    auto inputs = get_cell_inputs(pCell);
    auto outputs = static_cast<NeuralNetworkIntracellular*>(pCell->phenotype.intracellular)->forward(inputs);
    
    // Save inputs to custom_data
    for (size_t i = 0; i < inputs.size(); ++i) {
        pCell->custom_data["input_" + std::to_string(i)] = inputs[i];
//        std::cout << "Cell " << pCell->ID << " updated input_" << i << " to: " << inputs[i] << std::endl;

    }
    
    // Save outputs to custom_data
    for (size_t i = 0; i < outputs.size(); ++i) {
        pCell->custom_data["output_" + std::to_string(i)] = outputs[i];
//        std::cout << "Cell " << pCell->ID << " updated output_" << i << " to: " << outputs[i] << std::endl;

    }
    
////     Print inputs and outputs for debugging
//    std::cout << "Cell " << pCell->ID << " Inputs: ";
//    for (double input : inputs) {
//        std::cout << input << " ";
//    }
//    std::cout << std::endl;
//    
//    std::cout << "Cell " << pCell->ID << " Outputs: ";
//    for (double output : outputs) {
//        std::cout << output << " ";
//    }
//    std::cout << std::endl;

    
    // Interpret outputs and update custom data:
    pCell->custom_data["divide"] = (outputs[0] > 0.5) ? 1.0 : 0.0;
//    pCell -> phenotype.motility.migration_speed = 2*outputs[1];
    pCell->custom_data["separate"] = (outputs[1] > 0.5) ? 1.0 : 0.0;
    

    
    // Division
    if (pCell->custom_data["divide"] == 1.0) {
        // Disable motility
        pCell->phenotype.motility.is_motile = false;
        // Evaluate metabolism based on inactivity
        inactive_prey_growth_and_metabolism(pCell, phenotype, dt);
        
        //Inactive metebolic rate

        // Check if the cell has reached the size threshold for division
        double initial_volume = pCell->custom_data["initial_volume"];
        
        if (pCell->phenotype.volume.total >= 2 * initial_volume)
        {
            // Remove attachments and divide
            pCell->remove_all_attached_cells();
            pCell->remove_all_spring_attachments();
            
            // Cast the intracellular model to NeuralNetworkIntracellular
            NeuralNetworkIntracellular* nn_model = static_cast<NeuralNetworkIntracellular*>(pCell->phenotype.intracellular);

            // Apply mutation to the neural network
            // mutate(double mutation_rate, double weight_mutation_prob, double step_size)
            nn_model->mutate(1.0, 0.2, 0.1);  // Adjust mutation parameters as needed

            // Flag for division
            pCell-> phenotype.flagged_for_division = true;

        }
    }else{
        pCell->phenotype.motility.is_motile = true;
        active_prey_growth_and_metabolism(pCell, phenotype,dt);
    }
    
    
    // Seperation
    if (pCell->custom_data["separate"] == 1.0){
        // Remove attachments
        pCell -> remove_all_attached_cells();
    }

    
    // Predator stress evaluation
//    prey_stress_function(pCell, phenotype, dt);

    // Death evaluation
    evaluate_death_conditions(pCell, phenotype, dt);
    
    return; }
 
void inactive_prey_growth_and_metabolism(Cell* pCell, Phenotype& phenotype, double dt)
{
    
    // Update Energy:
    // No food consumption or growth posible but lower metabolic rate
    static double metabolic_rate = pCell -> custom_data["active_metabolic_rate"];
    pCell->custom_data["energy"] /= (1.0 + dt*metabolic_rate);
    
}


void active_prey_growth_and_metabolism(Cell* pCell, Phenotype& phenotype, double dt)
{
    // Access index for food
    int food_index = microenvironment.find_density_index("food");
    // sample the microenvironment at the cell’s locaiton
    double food_nearby = pCell->nearest_density_vector()[food_index];
    // Define the threshold and calculate food_eaten
    double threshold = 10.0;  // Example threshold
    double food_eaten;
    if (food_nearby < threshold)
        food_eaten = food_nearby;  // Eat all available food if below threshold
    else
        food_eaten = 10.0;  // Eat a fixed amount otherwise
    
    // Update Energy:
    // Cell gains energy equivalent to food consumption
    pCell->custom_data["energy"] += food_eaten; // Cell gains energy equivalent to food consumption
    // Update energy based on metabolism
    static double metabolic_rate = pCell -> custom_data["active_metabolic_rate"];
    pCell->custom_data["energy"] /= (1.0 + dt*metabolic_rate);
    
    // Update the microenvironment:
    // reduce food at the cell's location (it "ate" it)
    double new_food_level = std::max(0.0, food_nearby - food_eaten); // Ensure not negative
    int voxel_index = pCell->get_current_voxel_index();
    microenvironment(voxel_index)[food_index] = new_food_level;
    
    
    // ~~~~~~~~~~~~~~~~~~~~
    // Update Volume:
    // ~~~~~~~~~~~~~~~~~~~~
    
    // Power Law Model
    // Constants for the power law model
    const double b = 3/4 ; // Define the exponent b
    const double k_base = 10;  // Base value for k
    double k = k_base * food_eaten;

    // Get the current volume
    double current_volume = pCell->phenotype.volume.total;
    // Calculate the rate of change of volume
    double dVdt = k * pow(current_volume, b);
    // Update the volume based on the rate of change and the time step
    double new_volume = current_volume + dVdt * dt;
    // Set the new volume
    pCell->set_total_volume(new_volume);
}

void prey_stress_function(Cell* pCell, Phenotype& phenotype, double dt) {
    static Cell_Definition* pPreyDef = find_cell_definition("Cell");
    static Cell_Definition* pPredDef = find_cell_definition("Predator");

    
    // Variable to track if any predators are nearby
    bool predator_nearby = false;

    // Check for predators and increment danger level if a predator is nearby
    for (int n = 0; n < pCell->state.neighbors.size(); n++) {
        
        Cell* nbr = pCell->state.neighbors[n];
        if (nbr->type == pPredDef->type) {
            pCell->custom_data["danger"] += 1;
            predator_nearby = true;
//            std::cout << "Cell ID (" << pCell->ID << ") detects predator. " << std::endl;
        }
    }
    // If no predators are detected, decrease the danger level
     if (!predator_nearby) {
         pCell->custom_data["danger"] -= 0.3;
     }
    
    // Ensure danger does not go negative
    if (pCell->custom_data["danger"] < 0) {
        pCell->custom_data["danger"] = 0;
    }

}



void custom_division_function( Cell* pCell1, Cell* pCell2 )
{
    // Record division time (roughly)
    pCell1->custom_data["last_division_time"] = PhysiCell_globals.current_time;
    pCell2->custom_data["last_division_time"] = PhysiCell_globals.current_time;


    attach_cells(pCell1, pCell2);
    
    pCell1->custom_data["edible"] = 0.0;
    pCell2->custom_data["edible"] = 0.0;
    

    double initial_volume_1 = pCell1->custom_data["initial_volume"];
    double current_volume_1 = pCell1 ->phenotype.volume.total;
    
    double initial_volume_2 = pCell2->custom_data["initial_volume"];
    double current_volume_2 = pCell2 ->phenotype.volume.total;

    if ( current_volume_1 > 2*initial_volume_1 )
    {
        pCell1->divide();
        pCell2->divide();
    }

    return;
}

    
void contact_function( Cell* pMe, Phenotype& phenoMe , Cell* pOther, Phenotype& phenoOther , double dt )
{ return; } 

    
//// *********************************************************************************************
//// *********************************************************************************************
//// PREDATOR FUNCTONS
//// *********************************************************************************************
//// *********************************************************************************************
//
void pred_phenotype_function( Cell* pCell, Phenotype& phenotype, double dt ){
    
    
    // Hunt and eat prey cells (and increase energy)
    pred_hunt_function(pCell, phenotype, dt);
    
    // Update energy based on metabolism
    static double metabolic_rate = pCell -> custom_data["metabolic_rate"];
    pCell->custom_data["energy"] /= (1.0 + dt*metabolic_rate);
    
    
//    // Death based on starvation
//    if( pCell->custom_data["energy"] < 1.0 )
//    {
//        pCell->lyse_cell();
////        std::cout << "Predator Cell " << pCell->ID << ": Is dying" << std::endl;
//        
//    }
}
//
void pred_hunt_function( Cell* pCell, Phenotype& phenotype, double dt ){
    static Cell_Definition* pPreyDef = find_cell_definition( "Cell" );
    static Cell_Definition* pPredDef = find_cell_definition( "Predator" );
    
    bool is_hungry = (pCell->custom_data["energy"] < 800);
    double ingest_probability = pCell -> custom_data["prey_ingest_prob"];
    
    for( int n=0; n < pCell->state.neighbors.size(); n++ )
    {
        Cell* nbr = pCell->state.neighbors[n];
        
        if( nbr->type == pPreyDef->type && nbr->custom_data["edible"] == true && is_hungry ){
            
            double random_number = UniformRandom(); // Generate a random number between 0 and 1
            
            if (random_number < ingest_probability) {
                pCell -> ingest_cell(nbr);
                pCell -> custom_data["energy"] +=10;
                
                is_hungry = (pCell->custom_data["energy"] < 800);  // Re-evaluate hunger in case the condition changes due to energy increase
            }
        }
    }
}


