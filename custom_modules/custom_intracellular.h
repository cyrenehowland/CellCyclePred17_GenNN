//
//  custom_intracellular.h
//  
//
//  Created by Cyrene Howland on 2024-07-10.
//

#ifndef CUSTOM_INTRACELLULAR_h_
#define CUSTOM_INTRACELLULAR_h_

#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <iomanip>   // for setw

#include "/Users/cyrenehowland/PhysiCell/core/PhysiCell.h"
#include "/Users/cyrenehowland/PhysiCell/BioFVM/BioFVM.h"
#include "/Users/cyrenehowland/PhysiCell/core/PhysiCell_phenotype.h"
#include "/Users/cyrenehowland/PhysiCell/core/PhysiCell_cell.h"
#include "/Users/cyrenehowland/PhysiCell/modules/PhysiCell_pugixml.h"


// Define new class NeuralNetworkIntracellular inherited from PhysiCell::Intracellular
class NeuralNetworkIntracellular : public PhysiCell::Intracellular {
public:
    // declare 2D vector to store weights between input and hidden layer
    std::vector<std::vector<double>> weights_input_hidden;
    // declare 2D vector to store weights between hidden and output layer
    std::vector<std::vector<double>> weights_hidden_output;
    // declare a vector to store the activations of the hidden layer
    std::vector<double> hidden_layer;
    // declare a vector to store the activations of the output layer
    std::vector<double> output_layer;
    bool is_initialized;
    
    // Constructor, initialize neural network with the given sizes for input, hidden, and output layers
    NeuralNetworkIntracellular(int input_size, int hidden_size, int output_size);
    // initialize weights of NN
    void initialize_random(int seed);
    // perform a forward pass through NN and return output activations
    std::vector<double> forward(std::vector<double>& inputs);
    // mutate weights
    void mutate(double mutation_rate, double weight_mutation_prob, double step_size);
    
    // Intracellular interface methods:
    // initialize the intracellular model from an XML node
    void initialize_intracellular_from_pugixml(pugi::xml_node& node) override;
    // start the intracellular model
    void start() override;
    // check if the intracellular model needs updating
    bool need_update() override;
    // update the intracellular model
    void update() override;
    // update the cell's phenotype based on the intracellular model
    void update(PhysiCell::Cell* cell, PhysiCell::Phenotype& phenotype, double dt) override;
    // inherit properties from the parent cell
    void inherit(PhysiCell::Cell* cell) override;
    // get a parameter value by name
    double get_parameter_value(std::string name) override;
    // set a parameter value by name
    void set_parameter_value(std::string name, double value) override;
    // get the state of the intracellular model as a string.
    std::string get_state() override;
    // display the state of the intracellular model.
    void display(std::ostream& os) override;
    //  clone the intracellular model
    Intracellular* clone() override;
    
    bool has_variable(std::string name) override;
    bool get_boolean_variable_value(std::string name) override;
    void set_boolean_variable_value(std::string name, bool value) override;
    void print_current_nodes() override;
    int update_phenotype_parameters(PhysiCell::Phenotype& phenotype) override;
    int validate_PhysiCell_tokens(PhysiCell::Phenotype& phenotype) override;
    int validate_SBML_species() override;
    int create_custom_data_for_SBML(PhysiCell::Phenotype& phenotype) override;
    
    void print_weights() const;
    
};
// initialize the neural network for a given cell.
void initialize_cell_network(PhysiCell::Cell* cell, int input_size, int hidden_size, int output_size);
// get the inputs for the neural network from the cell
std::vector<double> get_cell_inputs(PhysiCell::Cell* cell);


#endif // CUSTOM_INTRACELLULAR_h_
