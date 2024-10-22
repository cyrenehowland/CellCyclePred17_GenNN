
// contains decleration of NeuralNetworkIntracellular class and methods
#include "./custom_intracellular.h"
// PhysiCell_phenotype class, which defines the structure and behavior of cells in the PhysiCell framework
#include "../core/PhysiCell_phenotype.h"
#include <cmath>
#include <cstdlib>
#include <ctime>


// Constructor for NeuralNetworkIntracellular class,
// Initialize a new instance of the class with sizes for the input layer, hidden layer, and output layer
NeuralNetworkIntracellular::NeuralNetworkIntracellular(int input_size, int hidden_size, int output_size)
    
    // Initialize weights with the appropriate sizes:
    : is_initialized(false) {
    
    // Initialize weight matrix between input and hidden layer w/ dims [input_size][hidden_size] filled w/ 0s (default)
    weights_input_hidden = std::vector<std::vector<double>>(input_size, std::vector<double>(hidden_size));
    // Initialize weight matrix between hidden layer and output layer w/ default values
    weights_hidden_output = std::vector<std::vector<double>>(hidden_size, std::vector<double>(output_size));
    // Initialize hidden layer and output layer
    hidden_layer = std::vector<double>(hidden_size);
    output_layer = std::vector<double>(output_size);
  
}

// Initialize random weights for neural network
void NeuralNetworkIntracellular::initialize_random(int seed) {
    srand(seed); // ensures different random values everytime the program runs
    
    // Iterate over each layer in the input-to-hidden weight matrix
    for (auto& layer : weights_input_hidden) {
        // Iterate over each weight in the current layer
        for (auto& weight : layer) {
            // Set each weight to a random value between -1 and 1
            weight = static_cast<double>(rand()) / RAND_MAX * 2.0 - 1.0;
        }
    }
    // Iterate over each layer in the hidden-to-output weight matrix
    for (auto& layer : weights_hidden_output) {
        for (auto& weight : layer) {
            // Set each weight to a random value between -1 and 1
            weight = static_cast<double>(rand()) / RAND_MAX * 2.0 - 1.0;
        }
    }
    is_initialized = true;
}

// Forward Pass: Computes outputs besed on inputs
std::vector<double> NeuralNetworkIntracellular::forward(std::vector<double>& inputs) {
    
    // Check if the object itself is valid
    if (this == nullptr) {
        std::cerr << "Error: The NeuralNetworkIntracellular object is null." << std::endl;
        return std::vector<double>();
    }
    
    if (!is_initialized) {
        std::cerr << "Error: Neural network is not initialized." << std::endl;
        return std::vector<double>();
    }
    // Ensure weights and layers are properly sized
    if (weights_input_hidden.empty() || weights_hidden_output.empty() ||
        hidden_layer.size() != weights_input_hidden[0].size() ||
        output_layer.size() != weights_hidden_output[0].size()) {
        std::cerr << "Error: Neural network layers or weights are not properly initialized or sized." << std::endl;
        return std::vector<double>(); // Return an empty vector in case of error
    }
    
    // Hidden layer calculation:
    // iterate over each neuron in the hidden layer
    for (size_t i = 0; i < hidden_layer.size(); ++i) {
        // initialize the hidden layer neuron value to zero
        hidden_layer[i] = 0.0;
        // iterate over each input
        for (size_t j = 0; j < inputs.size(); ++j) {
            // compute weighted sum of inputs for hidden layer neuron
            hidden_layer[i] += inputs[j] * weights_input_hidden[j][i];
        }
        hidden_layer[i] = std::max(0.0, hidden_layer[i]); // ReLU activation
    }

    // Output layer calculation:
    // iterate over each neuron in output layer
    for (size_t i = 0; i < output_layer.size(); ++i) {
        output_layer[i] = 0.0; // initialize output layer neuron to zero
        
        // iterate over each hidden layer neuron
        for (size_t j = 0; j < hidden_layer.size(); ++j) {
            // Check if indices are within bounds
            if (j >= weights_hidden_output.size() || i >= weights_hidden_output[j].size()) {
                std::cerr << "Error: Accessing weights_hidden_output out of bounds." << std::endl;
                return std::vector<double>();
            }
            // computes weighted sum of hidden layer outputs for the output layer neuron
            output_layer[i] += hidden_layer[j] * weights_hidden_output[j][i];
        }
        output_layer[i] = 1.0 / (1.0 + exp(-output_layer[i])); // Sigmoid activation
    }
    return output_layer;
}


// Print weights
void NeuralNetworkIntracellular::print_weights() const {
    std::cout << "Input to Hidden Weights:" << std::endl;
    for (const auto& row : weights_input_hidden) {
        for (double weight : row) {
            std::cout << weight << " ";
        }
        std::cout << std::endl;
    }

    std::cout << "Hidden to Output Weights:" << std::endl;
    for (const auto& row : weights_hidden_output) {
        for (double weight : row) {
            std::cout << weight << " ";
        }
        std::cout << std::endl;
    }
}


// Implement the required Intracellular interface methods here...
void NeuralNetworkIntracellular::initialize_intracellular_from_pugixml(pugi::xml_node& node) {} // Placeholder
void NeuralNetworkIntracellular::start() {} // Placeholder
bool NeuralNetworkIntracellular::need_update() { return true; } // Indicates that the intracellular model needs updating
void NeuralNetworkIntracellular::update() {} // Placeholder

// This method updates the state of the cell based on the neural network's output
void NeuralNetworkIntracellular::update(PhysiCell::Cell* cell, PhysiCell::Phenotype& phenotype, double dt) {
    // retrive inputs for neural network from the cell
    auto inputs = get_cell_inputs(cell);
    // compute output of neural network based on inputs
    auto outputs = forward(inputs);
//    // checks if the first output exceeds threshold
//    if (outputs[0] > 0.5) {
//        cell->custom_data["separate"] = 1.0;
//    } else {
//        cell->custom_data["separate"] = 0.0;
//    }
//    // checks if second output exceeds threshold
//    if (outputs[1] > 0.5) {
//        cell->custom_data["divide"] = 1.0;
//    } else {
//        cell->custom_data["divide"] = 0.0;
//    }
}



void NeuralNetworkIntracellular::mutate(double mutation_rate, double weight_mutation_prob, double step_size) {
    // Print the weights for verification
//   std::cout << "Mutating weights of parent cell.... " << std::endl;
//    std::cout << "Weights of parent cell before mutation " << std::endl;
//    print_weights();
    // Mutate input to hidden weights
    for (auto& layer : weights_input_hidden) {
        for (auto& weight : layer) {
            if (static_cast<double>(rand()) / RAND_MAX < mutation_rate) {
                if (static_cast<double>(rand()) / RAND_MAX < weight_mutation_prob) {
                    double original_weight = weight;
                    if (rand() % 2 == 0) {
                        weight += step_size;  // 50% chance to add step_size
                    } else {
                        weight -= step_size;  // 50% chance to subtract step_size
                    }

                }
            }
        }
    }

    // Mutate hidden to output weights
    for (auto& layer : weights_hidden_output) {
        for (auto& weight : layer) {
            if (static_cast<double>(rand()) / RAND_MAX < mutation_rate) {
                if (static_cast<double>(rand()) / RAND_MAX < weight_mutation_prob) {
                    double original_weight = weight;
                    if (rand() % 2 == 0) {
                        weight += step_size;  // 50% chance to add step_size
                    } else {
                        weight -= step_size;  // 50% chance to subtract step_size
                    }
                }
            }
        }
    }
    
    // Print the weights after mutation
//       std::cout << "Weights after mutation:" << std::endl;
//       print_weights();
}






// This method is called when a cell inherits properties from its parent.
void NeuralNetworkIntracellular::inherit(PhysiCell::Cell* cell) {
        auto* parent_nn = static_cast<NeuralNetworkIntracellular*>(cell->phenotype.intracellular);
        this->weights_input_hidden = parent_nn->weights_input_hidden;
        this->weights_hidden_output = parent_nn->weights_hidden_output;
        this->hidden_layer = parent_nn->hidden_layer;
        this->output_layer = parent_nn->output_layer;
        
        
        // Print the weights for verification
//       std::cout << "Inheriting weights from parent cell " << cell->ID << " to new cell." << std::endl;
//       print_weights();
    
//        // Apply mutation to the neural network
//       double mutation_rate = 1.0; // Probability of mutating a specific layer
//       double weight_mutation_prob = 0.2; // Example probability of mutating a specific weight
//       double step_size = 0.1; // Example step size
//       mutate(mutation_rate, weight_mutation_prob, step_size);

       // Print the weights after mutation
//       std::cout << "Weights after mutation:" << std::endl;
//       print_weights();
    }


// This method is used to get the value of a parameter by name.
// Here, it always returns 0.0, indicating no parameters are currently being managed by name.
double NeuralNetworkIntracellular::get_parameter_value(std::string name) { return 0.0; }
// This method is used to set the value of a parameter by name.
// It is empty, indicating no parameters are currently being managed by name.
void NeuralNetworkIntracellular::set_parameter_value(std::string name, double value) {}
// This method returns the state of the intracellular model as a string.
// It returns an empty string, indicating no state information is provided.
std::string NeuralNetworkIntracellular::get_state() { return ""; }
// This method is used to display the state of the intracellular model.
// It is empty, indicating no display behavior is defined.
void NeuralNetworkIntracellular::display(std::ostream& os) {}

// This method creates and returns a clone of the intracellular model.
// It uses the copy constructor to create a new instance of NeuralNetworkIntracellular.
PhysiCell::Intracellular* NeuralNetworkIntracellular::clone() {
    // Print the weights for verification
//   std::cout << "Cloning weights from parent cell " << cell->ID << " to new cell." << std::endl;
    return new NeuralNetworkIntracellular(*this);
}

// This method checks if a variable with the given name exists in the intracellular model.
// It always returns false, indicating no variables are managed by name.
bool NeuralNetworkIntracellular::has_variable(std::string name) { return false; }
// This method retrieves the value of a boolean variable by name.
// It always returns false, indicating no boolean variables are managed.
bool NeuralNetworkIntracellular::get_boolean_variable_value(std::string name) { return false; }
// This method sets the value of a boolean variable by name.
// It is empty, indicating no boolean variables are managed.
void NeuralNetworkIntracellular::set_boolean_variable_value(std::string name, bool value) {}
// This method prints the current state of the network nodes.
// It is empty, indicating no such printing behavior is defined.
void NeuralNetworkIntracellular::print_current_nodes() {}
// This method updates the phenotype parameters based on the intracellular model.
// It always returns 0, indicating no parameters are updated.
int NeuralNetworkIntracellular::update_phenotype_parameters(PhysiCell::Phenotype& phenotype) { return 0; }
// This method validates the PhysiCell tokens against the intracellular model.
// It always returns 0, indicating no validation is performed.
int NeuralNetworkIntracellular::validate_PhysiCell_tokens(PhysiCell::Phenotype& phenotype) { return 0; }

// This method validates the SBML species against the intracellular model.
// It always returns 0, indicating no validation is performed.
int NeuralNetworkIntracellular::validate_SBML_species() { return 0; }
// This method creates custom data for SBML based on the intracellular model.
// It always returns 0, indicating no custom data is created.
int NeuralNetworkIntracellular::create_custom_data_for_SBML(PhysiCell::Phenotype& phenotype) { return 0; }
