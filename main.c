/*****************************************************************
 *  Kuramoto Model Simulation                                   *
 *  This program simulates synchronization in a network using a *
 *  second-order Kuramoto model with community detection.       *
 *****************************************************************/

#include <time.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <sys/time.h>

#include "src/utils/sampler.h"
#include "src/utils/normalizer.h"
#include "src/utils/solvers/rk4.h"
#include "src/utils/solvers/rkf45.h"
#include "src/utils/type_handlers.h"
#include "src/utils/solvers/gbs_integrator.h"
#include "src/utils/random/random_generators.h"

/*****************************************************************
*      Definitions                                              *
*****************************************************************/

/* ---- Constants ---- */
#define PI2 (2.0 * M_PI)  

/* ---- Randomization Parameters ---- */
int     seed1, seed3;           
unsigned long unique_id;        

/* ---- Network Parameters ---- */
int     num_nodes;        
int     num_edges;        
int     initial_edge_cut; 
int     failure_count;    
int     num_cut_nodes;    
char    network_type[250];
int     random_self_freq;
int     force_freq_from_file;

/* ---- Synchronization Model Parameters ---- */
double  initial_phase_randomness;   
double  coupling_0, coupling_target;
double  coupling_current;     
double  damping_coefficient;  
double  failure_threshold;    
double  precision_epsilon;   

/* HVDC parameters */
double  default_power;        

/* ---- Time & Iteration Control ---- */
int     num_stages;           
int     max_it_thermal;   
int     max_it_cascade;   
double  log_base;         
double  delta;            

/* ---- Time Evolution Tracking ---- */
int     *time_thermal, *time_cascade;
int     max_time_thermal, max_time_cascade;
double  *lambda_thermal;

/* ---- Edge Structure ---- */
double  total_coupling_weight; 
typedef struct {
    int nodes[2];       
    int edge_type;      
    double weight;   
} Edge;

Edge** adj_list;   

int symmetric_network;

/* ---- Node Structure ---- */
typedef struct {
  int node_id;       
  int degree;        
  int community;     
  double frequency;  
  int* neighbors;    
} Node;

Node* nodes;

/* ---- Network Representation ---- */
double  *state_variables;
double  *derivatives;


/* ---- Output Files ---- */
char    io_dir[250];
char    edges_filename[550];
char    nodes_filename[550];

/* ---- Scene Identifier ---- */
int     scene;       

/* ---- Global History Tracking Arrays ---- */
double *ord_param_hist_th;
double *freq_spread_hist_th;
double *univ_ord_param_hist_th;

double *local_ord_param;
double *local_freq_spread;
double *local_phase;

double *ord_param_hist_cascade;
double *freq_spread_hist_cascade;
double *univ_ord_param_hist_cascade;

int *failure_counts;
double *failure_times;
int **failure_nodes;  
double **last10Frequencies; 

/* ---- Community Detection ---- */
int do_community; 
int num_communities;

/* ---- Community Tracking ---- */
double **comm_ord_param_hist_therm;
double **comm_univ_ord_param_hist_therm;
double **comm_freq_spread_hist_therm;

double **comm_ord_param_hist_cascade;
double **comm_univ_ord_param_hist_cascade;
double **comm_freq_spread_hist_cascade;

double *comm_sum_ord_param;
double *comm_sum_univ_ord_param;
double *comm_sum_freq_sum;
double *comm_sum_freq_square_sum;
double *comm_weights; 
int *comm_count;

/* ---- Printing/Info ---- */
int verbose = 0; 
int quiet_run = 0;

/* ---- Placeholder functions ---- */
void print_help(void);
void main_program(void);
void generate_network(void);
void initialize_system(void);
void pre_init_alloc(void);
void post_init_alloc(void);
void step(void);
void write_output(void);

/**
 * @brief Main function to initialize parameters, parse arguments, and start the simulation.
 * 
 * @param argc Number of command-line arguments.
 * @param argv Array of command-line arguments.
 * @return int Returns 0 upon successful execution.
 */
int main(int argc, char **argv)
{
    int i;
    long i1;
    struct timeval tv;
    int show_help = 0;

    /* ---- Default Parameter Initialization ---- */
    seed1 = 1976;
    num_edges = 0;
    num_nodes = 0;
    delta = 0.01;
    symmetric_network = 1;
    random_self_freq = 1;
    force_freq_from_file = 0;
    precision_epsilon = 1.0e-12;
    max_it_thermal = 1000;
    max_it_cascade = 1000;
    default_power = 0.5;
    coupling_0 = 0.01;
    coupling_target = 0.5;
    damping_coefficient = 0.1;
    initial_phase_randomness = 0.0;
    log_base = 1.08;
    failure_threshold = 0.7;
    initial_edge_cut = -1;
    num_cut_nodes = 1;
    num_stages = 10;  
    do_community = 0; 
    strcpy(io_dir, "hu387");
    strcpy(network_type, "base");

    /* ---- Read Command-line Arguments ---- */
    for (i = 1; i < argc; i++)
    {
        if (strcmp(argv[i], "--h") == 0) {
            print_help();
            exit(0);
        }
        
        char *value = strchr(argv[i], '=');
        if (value != NULL) {
            *value = '\0';
            value++; 
        } else {
            printf("Error: Invalid argument format. Use -option=value or --option=value\n");
            exit(1);
        }

        if (strcmp(argv[i], "-seed") == 0) seed1 = atoi(value);
        else if (strcmp(argv[i], "-init_rand") == 0) initial_phase_randomness = atof(value);
        else if (strcmp(argv[i], "-k") == 0) coupling_target = atof(value);
        else if (strcmp(argv[i], "-k0") == 0) coupling_0 = atof(value);
        else if (strcmp(argv[i], "-log_base") == 0) log_base = atof(value);
        else if (strcmp(argv[i], "-def_pow") == 0) default_power = atof(value);
        else if (strcmp(argv[i], "-it_th") == 0) max_it_thermal = atoi(value);
        else if (strcmp(argv[i], "-it_lc") == 0) max_it_cascade = atoi(value);
        else if (strcmp(argv[i], "-alpha") == 0) damping_coefficient = atof(value);
        else if (strcmp(argv[i], "-io_dir") == 0) strcpy(io_dir, value);
        else if (strcmp(argv[i], "-delta") == 0) delta = atof(value);
        else if (strcmp(argv[i], "-stages") == 0) num_stages = atoi(value);
        else if (strcmp(argv[i], "-sym") == 0) symmetric_network = atoi(value);
        else if (strcmp(argv[i], "-rsf") == 0) random_self_freq = atoi(value);
        else if (strcmp(argv[i], "-ff") == 0) force_freq_from_file = atoi(value);
        else if (strcmp(argv[i], "-ic") == 0) initial_edge_cut = atoi(value);
        else if (strcmp(argv[i], "-vb") == 0) verbose = atoi(value);
        else if (strcmp(argv[i], "-quiet") == 0) quiet_run = atoi(value);
        else if (strcmp(argv[i], "-thr") == 0) failure_threshold = atof(value);
        else if (strcmp(argv[i], "-scene") == 0) scene = atoi(value);
        else if (strcmp(argv[i], "-id") == 0) unique_id = atoi(value);
        else if (strcmp(argv[i], "-nt") == 0) strcpy(network_type, value);
        else if (strcmp(argv[i], "-comms") == 0) do_community = atoi(value);
    }
    int max_time_thermal_tmp = (int)(log((double)max_it_thermal)/log(log_base))+1;
    max_time_thermal = (max_time_thermal_tmp < max_it_thermal) ? max_time_thermal_tmp : max_it_thermal;
    int max_time_cascade_tmp = (int)(log((double)max_it_cascade)/log(log_base))+1;
    max_time_cascade = (max_time_cascade_tmp < max_it_cascade) ? max_time_cascade_tmp : max_it_cascade;
        
    i1 = gettimeofday(&tv, NULL);
    // if (!unique_id) unique_id = tv.tv_usec;
    
    unsigned long uunique_id = tv.tv_usec;
    seed1 = uunique_id;
    seed3 = -time(NULL);

    if (!verbose) printf("\n# Allocating memory pre init...\n");
    pre_init_alloc();
    
    if (!verbose) printf("# Initializing system...\n");
    initialize_system();

    if (!quiet_run) {
        printf("\n=== Simulation Parameters ===\n");
        printf("# System size (N)       : %d\n", num_nodes);
        printf("# Number of edges (E)   : %d\n", num_edges);
        printf("# Thermalization steps  : %d\n", max_it_thermal);
        printf("# Cascade steps         : %d\n", max_it_cascade);
        printf("# Damping coefficient   : %lf\n", damping_coefficient);
        printf("# Target Coupling       : %lf\n", coupling_target);
        printf("# Initial Randomness    : %lf\n", initial_phase_randomness);
        printf("# Failure threshold     : %lf\n", failure_threshold);
        printf("# Cut edges             : %d\n", initial_edge_cut);
        printf("# Symmetric network     : %d\n", symmetric_network);
        printf("# Random self-freq      : %d\n", random_self_freq);
        printf("# Force freq from file  : %d\n", force_freq_from_file);
        printf("# Network type          : %s\n", network_type);
        printf("# Community analysis    : %d\n", do_community);
    }
    
    if (!verbose) printf("\n# Allocating memory post init...\n");
    post_init_alloc();

    if (!verbose) printf("\n# Generating network...\n");
    generate_network();

    if (!verbose) printf("\n# Starting simulation...\n");
    step();
    
    if (!verbose) printf("\n# Writing output...\n");
    write_output();
    
    if (!verbose) printf("\n# Simulation complete.\n");

    return 0;
}


/**
 * @brief Displays help information about command-line arguments.
 * 
 * This function prints a list of available command-line options 
 * along with their descriptions.
 */
void print_help(void)
{
    printf("\nKuramoto Model Simulation - Help Menu\n");
    printf("=====================================\n");
    printf("Usage: main [OPTIONS]\n\n");

    printf("---- Network Configuration ----\n");
    printf("  -nt[type]    : Network topology type (default: 'base')\n");
    printf("  -scene[num]  : Scene ID\n");
    printf("  -comms[bool]    : Do community analysis\n");

    printf("\n---- Model Parameters ----\n");
    printf("  -seed[num]       : Set random seed\n");
    printf("  -init_rand[num]  : Initial phase randomness\n");
    printf("  -k[num]          : Target coupling value\n");
    printf("  -k0[num]         : Initial coupling value\n");
    printf("  -alpha[num]      : Damping coefficient (alpha)\n");
    printf("  -thr[num]        : Failure threshold\n");
    printf("  -def_pow[num]    : Default power in/out flow at HVDC/DC lines\n");
    printf("  -log_base[num]   : Logarithmic base for output time intervals\n");
    printf("  -delta[num]      : Integration step delta\n");

    printf("\n---- Simulation Control ----\n");
    printf("  -it_th[num]  : Max thermalization steps\n");
    printf("  -it_lc[num]  : Max cascade steps\n");
    printf("  -stages[num] : Number of stages for thermalization\n");
    printf("  -sym[bool]  : Symmetric network edges\n");
    printf("  -rsf[bool]   : Random self-frequency assignment (0-centered Gaussian)\n");
    printf("  -ic[num]     : Initial edge ID to cut\n");
    printf("  -id[num]     : Unique ID\n");

    printf("\n---- Output Options ----\n");
    printf("  -io_dir[file]  : I/O directory name\n");
    
    printf("\n---- Miscellaneous ----\n");
    printf("  -vb[num]     : Verbose output\n");
    printf("  -quiet[num]  : Suppress output\n");
    printf("  --h          : Display this help menu\n");

    printf("\nExample:\n");
    printf("  ./main -seed=1234 -k=0.5 -alpha=0.1 -it_th=1000 -it_lc=1000 -io_dir=results -nt=base -scene=1\n\n");
}

/***************************************************************************
 * @brief Initializes the system parameters and simulation variables.
 *
 * This function:
 * - Reads node frequencies and community assignments from a file.
 * - Initializes order parameter tracking arrays.
 * - Sets up logarithmic binning for time intervals.
 ***************************************************************************/

void initialize_system() {
    double base_frequency, weight;
    int t, i, row, col, community_id, edge_type, force_freq_bool;
    force_freq_bool = 0;
    
    FILE *file;

    sprintf(edges_filename, "./networks/%s/%s/edges%d.data", io_dir, network_type, scene);

    file = fopen(edges_filename, "r");
    if (!file) {
        fprintf(stderr, "Error: Could not open edges file %s\n", edges_filename);
        exit(EXIT_FAILURE);
    }

    while (fscanf(file, "%d %d %lf %d\n", &row, &col, &weight, &edge_type) == 4) {
        num_edges++;
        if (row > num_nodes) num_nodes = row;
        if (col > num_nodes) num_nodes = col;
    }
    fclose(file);

    if (initial_edge_cut == -1) {
        srand(time(NULL));
        initial_edge_cut = rand() % num_edges;
    }

    sprintf(nodes_filename, "./networks/%s/%s/nodes%d.data", io_dir, network_type, scene);
    file = fopen(nodes_filename, "r");

    if (!file) {
        printf("Error: Unable to open node data file: %s\n", nodes_filename);
        exit(EXIT_FAILURE);
    }

    nodes = (Node*)malloc((num_nodes + 1) * sizeof(Node));
    if (!nodes) {
        printf("Error: Memory allocation failed for nodes.\n");
        exit(EXIT_FAILURE);
    }

    for (i = 1; i <= num_nodes; i++) {
        community_id = -1;

        if (force_freq_from_file) {
            char *format = "%d %lf %d %d\n";
            if (fscanf(file, format, &row, &base_frequency, &community_id, &force_freq_bool) < 4) {
                printf("Error: Incorrect file format or insufficient data in %s.\n", nodes_filename);
                exit(EXIT_FAILURE);
            }
        } else {
            char *format = "%d %lf %d\n";
            if (fscanf(file, format, &row, &base_frequency, &community_id) < 3) {
                printf("Error: Incorrect file format or insufficient data in %s.\n", nodes_filename);
                exit(EXIT_FAILURE);
            }
        }
        
        if (community_id > num_communities) num_communities = community_id;
        nodes[i].node_id = row;
        nodes[i].degree = 0;  // Will be set in `generate_network()`
        nodes[i].community = community_id;
        if (random_self_freq) {
            if (force_freq_from_file && force_freq_bool) {
                nodes[i].frequency = base_frequency;
            } else {            
                nodes[i].frequency = gasdev(&seed1);
            }
        } else {
            if (force_freq_from_file && force_freq_bool) {
                nodes[i].frequency = base_frequency;
            } else {            
                nodes[i].frequency = base_frequency + gasdev(&seed1) * base_frequency * 0.05;
            }
        }
        nodes[i].neighbors = NULL;  // Will be allocated in `generate_network()`
    }
    fclose(file);
    
    printf("\n# Node 0: Id: %d, Frequency: %lf, Community: %d\n", nodes[1].node_id, nodes[1].frequency, nodes[1].community);
    printf("# Node 100: Id: %d, Frequency: %lf, Community: %d\n", nodes[100].node_id, nodes[100].frequency, nodes[100].community);

    if(do_community) {
        comm_count = (int*)calloc(num_communities + 1, sizeof(double));
        memset(comm_count,0,(num_communities+1)*sizeof(int));
        for (i=1; i<=num_nodes; i++) {
            comm_count[nodes[i].community]++;
        }
    }

    sample_logarithmically(max_it_thermal, max_time_thermal, log_base, 0, time_thermal);
    sample_logarithmically(max_it_cascade, max_time_cascade, log_base, 0, time_cascade);

    if (num_stages == 0) {
        for (t = 1; t <= max_it_thermal; t++) {
            lambda_thermal[t] = coupling_target;
        }
    } else {
        double log_factor = pow(max_it_thermal, 1.0 / num_stages);
        int stage_start = 1;
        double lambda_increment = (coupling_target - coupling_0) / (num_stages - 1);
        if (!verbose) printf("Target lambda: %lf, lambda increment %lf\n", coupling_target, lambda_increment);
        for (int s = 0; s < num_stages; s++) {
            int stage_end = (int)(pow(log_factor, s + 1));
            if (stage_end > max_it_thermal) stage_end = max_it_thermal;

            for (t = stage_start; t <= stage_end; t++) {
                lambda_thermal[t] = coupling_0 + lambda_increment * s;
            }
            stage_start = stage_end + 1;
        }

        /* Ensure the remaining time intervals are filled */
        for (t = stage_start; t <= max_it_thermal; t++) {
            lambda_thermal[t] = coupling_target;
        }
    }
}


/***************************************************************************
 * @brief Generates the network dynamically by reading edges from a file.
 *
 * This function:
 * - Reads the edge list from `edges.dat`.
 * - Dynamically determines `num_nodes` and `num_edges`.
 * - Constructs the adjacency list representation, storing edges separately.
 * - If `do_community` is enabled, calculates intra-community weights.
 ***************************************************************************/

void generate_network() {
    double edge_weight;
    int node1, node2, edge_type;
    
    FILE* file;

    sprintf(edges_filename, "./networks/%s/%s/edges%d.data", io_dir, network_type, scene);
    file = fopen(edges_filename, "r");

    if (!file) {
        printf("Error: Unable to open edge file: %s\n", edges_filename);
        exit(EXIT_FAILURE);
    }

    adj_list = (Edge**)malloc((num_edges + 1) * sizeof(Edge*));

    if (!nodes || !adj_list) {
        printf("Error: Memory allocation failed for network storage.\n");
        exit(EXIT_FAILURE);
    }

    file = fopen(edges_filename, "r");

    if (!file) {
        printf("Error: Unable to reopen edge file: %s\n", edges_filename);
        exit(EXIT_FAILURE);
    }

    total_coupling_weight = 0.0; 
    int edge_index = 0;

    while (fscanf(file, "%d %d %lf %d\n", &node1, &node2, &edge_weight, &edge_type) == 4) {
        adj_list[edge_index] = (Edge*)malloc(sizeof(Edge));
        if (!adj_list[edge_index]) {
            printf("Error: Memory allocation failed for edge %d.\n", edge_index);
            exit(EXIT_FAILURE);
        }
        
        adj_list[edge_index]->nodes[0] = node1;
        adj_list[edge_index]->nodes[1] = node2;
        adj_list[edge_index]->edge_type = edge_type;
        adj_list[edge_index]->weight = edge_weight;

        nodes[node1].degree++;
        nodes[node2].degree++;
        
        nodes[node1].neighbors = (int*)realloc(nodes[node1].neighbors, nodes[node1].degree * sizeof(int));
        nodes[node2].neighbors = (int*)realloc(nodes[node2].neighbors, nodes[node2].degree * sizeof(int));
        
        if (!nodes[node1].neighbors || !nodes[node2].neighbors) {
            printf("Error: Memory reallocation failed for adjacency list.\n");
            exit(EXIT_FAILURE);
        }
        
        nodes[node1].neighbors[nodes[node1].degree - 1] = node2;
        nodes[node2].neighbors[nodes[node2].degree - 1] = node1;

        total_coupling_weight += edge_weight;
        edge_index++;
    }
    fclose(file);
    total_coupling_weight = 1.0 / total_coupling_weight;

    if (do_community) {
        memset(comm_weights, 0, (num_communities + 1) * sizeof(double));
        for (int i = 0; i < num_edges; i++) {
            int source = adj_list[i]->nodes[0];
            int target = adj_list[i]->nodes[1];

            if (nodes[source].community == nodes[target].community) {
                comm_weights[nodes[source].community] += adj_list[i]->weight;
            }
        }

        for (int i = 1; i <= num_communities; i++) {
            if (comm_weights[i] > 0) {
                comm_weights[i] = 1.0 / comm_weights[i];
            }
        }
    }

    if (!quiet_run) printf("\n# Detected Nodes: %d, Edges: %d, Symmetrize %d, Total weight: %lf\n", num_nodes, num_edges, symmetric_network, total_coupling_weight);
}

/***************************************************************************
 * @brief Computes the derivatives for the second-order Kuramoto model.
 *
 * This function calculates the rate of change of phase and frequency for
 * each node using the adjacency list representation.
 *
 * @param x     Current time (not used in computation but kept for solvers).
 * @param y     Current state vector (phase and frequency values).
 * @param dydx  Output derivative vector.
 ***************************************************************************/

void compute_derivatives(double x, double y[], double dydx[]) {
    int i, edge_type;
    double interaction, power_1, power_2, freq_diff;

    for (i = 1; i <= num_nodes; i++) {
        dydx[i] = y[i + num_nodes]; // Frequency is the first derivative of phase
        dydx[i+num_nodes] = nodes[i].frequency - damping_coefficient * y[i + num_nodes];
    }

    for (i = 0; i < num_edges; i++) {
        double weight = adj_list[i]->weight;
        int node1 = adj_list[i]->nodes[0];
        int node2 = adj_list[i]->nodes[1];
        edge_type = adj_list[i]->edge_type;
        power_1 = 0;
        power_2 = 0;
        interaction = weight * sin(y[node2] - y[node1]);
        freq_diff = y[node2 + num_nodes] - y[node1 + num_nodes];

        if (edge_type == 1) {

            // printf("Edge type %d\n", edge_type);
            // printf("Frequency 1: %f, Frequency 2: %f interaction %f \n", dydx[node1], dydx[node2], interaction);
            /*
            if (y[node1 + num_nodes] > y[node2 + num_nodes]) {
                power_1 = - y[node1 + num_nodes] / (fabs(y[node1 + num_nodes]) + 1) * default_power;
                power_2 = y[node2 + num_nodes] / (fabs(y[node2 + num_nodes])  + 1) * default_power;
            } else {
                power_1 = y[node1 + num_nodes] / (fabs(y[node1 + num_nodes]) + 1) * default_power;
                power_2 = - y[node2 + num_nodes] / (fabs(y[node2 + num_nodes]) + 1) * default_power;
            }
            // printf("Node1: %d, Node2: %d, Power1: %lf, Power2: %lf\n", node1, node2, power_1, power_2);
            */

            dydx[node1 + num_nodes] +=  norm_tanh(freq_diff) * weight;
            dydx[node2 + num_nodes] += - norm_tanh(freq_diff) * weight;
        } else {
            dydx[node1 + num_nodes] += coupling_current * interaction;
            dydx[node2 + num_nodes] += -coupling_current * interaction * (double)symmetric_network;
        }
    }
}

/***************************************************************************
 * @brief Computes derivatives with line failure detection.
 *
 * This function computes the phase derivatives while checking for failures.
 * If an edge exceeds the failure threshold, it is removed from the network.
 *
 * @param x        Current time (not used in computation but kept for solvers).
 * @param y        Current state vector (phase and frequency values).
 * @param dydx     Output derivative vector.
 * @param lambda   Coupling strength.
 * @param time     Current simulation time (used for failure logging).
 ***************************************************************************/

void compute_derivatives_with_failures(double x, double y[], double dydx[], double lambda, int time) {
    int i, edge_type;
    double interaction, power_flow, power_1, power_2, freq_diff;

    for (i = 1; i <= num_nodes; i++) {
        dydx[i] = y[i + num_nodes]; // Frequency is the first derivative of phase
        dydx[i+num_nodes] = nodes[i].frequency - damping_coefficient * y[i + num_nodes];
    }

    for (i = 0; i < num_edges; i++) {
        double weight = adj_list[i]->weight;
        int node1 = adj_list[i]->nodes[0];
        int node2 = adj_list[i]->nodes[1];
        edge_type = adj_list[i]->edge_type;
        power_1 = 0;
        power_2 = 0;

        if (weight == 0) continue;

        power_flow = sin(y[node2] - y[node1]);

        if (fabs(power_flow) > failure_threshold && weight != 0.0) {
            if (!verbose) printf("# Line failure detected at time %d between nodes %d and %d\n", time, node1, node2);

            adj_list[i]->weight = 0;
            failure_count++;
            failure_times[failure_count] = time;
            failure_nodes[failure_count][0] = node1;
            failure_nodes[failure_count][1] = node2;
        }

        interaction = power_flow * weight;
        if (edge_type == 1) {
            /*
            if (y[node1 + num_nodes] > y[node2 + num_nodes]) {
                power_1 = - y[node1 + num_nodes] / nodes[node1].frequency * default_power;
                power_2 = y[node2 + num_nodes] / nodes[node2].frequency * default_power;
            } else {
                power_1 = y[node1 + num_nodes] / nodes[node1].frequency * default_power;
                power_2 = - y[node2 + num_nodes] / nodes[node2].frequency * default_power;
            }
            */
            dydx[node1 + num_nodes] += norm_tanh(norm_tanh(freq_diff)) * weight;
            dydx[node2 + num_nodes] += - norm_tanh(norm_tanh(freq_diff)) * weight;
        } else {
            dydx[node1 + num_nodes] += lambda * interaction;
            dydx[node2 + num_nodes] += -lambda * interaction * (double)symmetric_network;
        }
    }
}

/***************************************************************************
 * @brief Performs the main simulation steps using r8_rkf45 for integration.
 *
 * This function:
 * - Integrates the system over time using RKF45 adaptive integration.
 * - Tracks synchronization order parameters (global & community).
 * - Detects and processes failures in cascade mode.
 * - Stores simulation results in memory to be written later.
 ***************************************************************************/
void step() {
    int i, t;
    double complex z;
    int MAX_SAVE = 10;
    int stage_count = 0;
    int output_count = 1;
    int failure_count_previous = 0;
    double size_scale = 1.0 / (double) num_nodes;
    double x, x_out, order_param, prev_order_param,  univ_order_param;
    double freq_spread_square, freq_spread, freq_spread_sum, node_frequency;
    
    double relerr = 1e-7, abserr = 1e-7;
    int flag = 1;
    
    z = 0;
    x = 0.0;
    x_out = 0.0;
    freq_spread_sum = 0;
    prev_order_param = 0;
    univ_order_param = 0;
    freq_spread_square = 0;

    //TODO: put this into a reset_state() fn
    for (i = 1; i <= num_nodes; i++) {
        state_variables[i] = initial_phase_randomness * PI2 * ran_lagged_fibonacci(&seed1);
        node_frequency = nodes[i].frequency;
        state_variables[i + num_nodes] = node_frequency;
        z += cexp(I * state_variables[i]);

        freq_spread_sum += node_frequency;
        freq_spread_square += node_frequency * node_frequency;
    }

    for(i = 0; i < num_edges; i++){
        univ_order_param += adj_list[i]->weight * cos(state_variables[adj_list[i]->nodes[1]] - state_variables[adj_list[i]->nodes[0]]);
    }
    univ_order_param *= total_coupling_weight;
    order_param = cabs(z) * size_scale;
    prev_order_param = 0;
    freq_spread_sum *= size_scale;
    freq_spread = freq_spread_square * size_scale - freq_spread_sum * freq_spread_sum;
    
    if (!quiet_run) printf("# Initial in therm. R: %lf Omega: %lf R_uni %lf\n", order_param, freq_spread, univ_order_param);

    /* ============================= */
    /* ===== THERMALIZATION ======== */
    /* ============================= */

    if (!quiet_run) printf("\n# Stage %d, lambda %lf\n", stage_count+1, lambda_thermal[1]);
    for (t = 1; t <= max_it_thermal; t++) {
        x_out = x + delta;
        coupling_current = lambda_thermal[t];

        if (t > 1 && lambda_thermal[t] != lambda_thermal[t - 1]) {
            stage_count++;
            if (!quiet_run) printf("# Stage %d, lambda = %lf\n", stage_count+1, coupling_current);
        }

        compute_derivatives(x, state_variables, derivatives);
        flag = r8_rkf45(compute_derivatives, num_nodes * 2, state_variables, derivatives, &x, x_out, &relerr, abserr, flag);
        // gbs_integrator(compute_derivatives, num_nodes * 2, state_variables, derivatives, &x, x_out, &relerr, abserr, &flag);
        x += delta;

        if(t == time_thermal[output_count-1]){
            z = 0;
            freq_spread_sum = 0;
            univ_order_param = 0;
            freq_spread_square = 0;
            
            if(do_community){
                memset(comm_sum_ord_param, 0, (num_communities + 1) * sizeof(double));
                memset(comm_sum_univ_ord_param, 0, (num_communities + 1) * sizeof(double));
                memset(comm_sum_freq_square_sum, 0, (num_communities + 1) * sizeof(double));
                memset(comm_sum_freq_sum, 0, (num_communities + 1) * sizeof(double));
            
                for (i = 1; i <= num_nodes; i++) {
                    z += cexp(I * state_variables[i]);
                    node_frequency = state_variables[i + num_nodes];
                    freq_spread_square += node_frequency * node_frequency;
                    freq_spread_sum += node_frequency;
                    
                    int comm_idx = nodes[i].community;
                    comm_sum_ord_param[comm_idx] += cexp(I * state_variables[i]);
                    comm_sum_freq_square_sum[comm_idx] += node_frequency * node_frequency;
                    comm_sum_freq_sum[comm_idx] += node_frequency;
                    
                    last10Frequencies[i][output_count % MAX_SAVE] = state_variables[i + num_nodes];
                }

                for(i = 0; i < num_edges; i++){
                    double weight = adj_list[i]->weight;
                    int node1 = adj_list[i]->nodes[0];
                    int node2 = adj_list[i]->nodes[1];
                    univ_order_param += weight * cos(state_variables[node2] - state_variables[node1]);

                    if (nodes[node1].community == nodes[node2].community) {
                        comm_sum_univ_ord_param[nodes[node1].community] += weight * cos(state_variables[node2] - state_variables[node1]);
                    }
                }

                for(i = 1; i <= num_communities; i++){
                    comm_ord_param_hist_therm[output_count][i] = cabs(comm_sum_ord_param[i]) / (double)comm_count[i];
                    comm_univ_ord_param_hist_therm[output_count][i] = comm_sum_univ_ord_param[i] * comm_weights[i];
                    comm_freq_spread_hist_therm[output_count][i] = comm_sum_freq_square_sum[i] / (double)comm_count[i] - (comm_sum_freq_sum[i] / comm_count[i]) * (comm_sum_freq_sum[i] / comm_count[i]);
                }
            } else {
                for (i = 1; i <= num_nodes; i++) {
                    z += cexp(I * state_variables[i]);
                    node_frequency = state_variables[i + num_nodes];
                    freq_spread_square += node_frequency * node_frequency;
                    freq_spread_sum += node_frequency;

                    last10Frequencies[i][output_count % MAX_SAVE] = state_variables[i + num_nodes];
                }

                for(i = 0; i < num_edges; i++){
                    double weight = adj_list[i]->weight;
                    int node1 = adj_list[i]->nodes[0];
                    int node2 = adj_list[i]->nodes[1];
                    univ_order_param += weight * cos(state_variables[node2] - state_variables[node1]);
                }
            }

            univ_order_param *= total_coupling_weight;

            order_param = cabs(z) * size_scale;
            freq_spread_sum = freq_spread_sum * size_scale;
            freq_spread = freq_spread_square * size_scale - freq_spread_sum * freq_spread_sum;

            ord_param_hist_th[output_count] = order_param;
            freq_spread_hist_th[output_count] = freq_spread;
            univ_ord_param_hist_th[output_count] = univ_order_param;
            
            if (!verbose) printf("%d %lf %lf %lf %lf\n",t, order_param, freq_spread, univ_order_param, coupling_current);

            if (t > 9000 && fabs(order_param - prev_order_param) < 0.0001 && coupling_current == coupling_target) {
                break;
            }
            output_count ++;
            prev_order_param = order_param;
        }

    }

    for (i = 1; i <= num_nodes; i++) {
        double complex order_param_z = 0;
        double loc_freq_sum_square = 0, loc_freq_sum = 0;
        double loc_freq_spread = 0;

        for (int k = 0; k < nodes[i].degree; k++) {
            int neighbor = nodes[i].neighbors[k];
            order_param_z += cexp(I * state_variables[neighbor]);
            loc_freq_sum += state_variables[neighbor + num_nodes];
            loc_freq_sum_square += state_variables[neighbor + num_nodes] * state_variables[neighbor + num_nodes];
        }

        local_ord_param[i] = cabs(order_param_z) / (double)nodes[i].degree;
        loc_freq_sum /= (double)nodes[i].degree;
        local_freq_spread[i] = loc_freq_sum_square / (double)nodes[i].degree - loc_freq_sum * loc_freq_sum;
        local_phase[i] = state_variables[i];
    }

    /* ============================== */
    /* ===== CASCADE FAILURE ======== */
    /* ============================== */

    z = 0;
    x = 0.0;
    output_count = 1;
    failure_count = 0;
    freq_spread_sum = 0;
    univ_order_param = 0;
    freq_spread_square = 0;
    coupling_current = coupling_target;

    if (!quiet_run) printf("\n Casc time %d %d %d %d\n", time_cascade[0], time_cascade[1], time_cascade[10], time_cascade[100]);

    for (i = 1; i <= num_nodes; i++) {
        node_frequency = state_variables[i + num_nodes];
        z += cexp(I * state_variables[i]);

        freq_spread_sum += node_frequency;
        freq_spread_square += node_frequency * node_frequency;
    }

    for(i = 0; i < num_edges; i++){
        univ_order_param += adj_list[i]->weight * cos(state_variables[adj_list[i]->nodes[1]] - state_variables[adj_list[i]->nodes[0]]);
    }

    prev_order_param = 0;
    freq_spread_sum *= size_scale;
    order_param = cabs(z) * size_scale;
    univ_order_param *= total_coupling_weight;
    freq_spread = freq_spread_square * size_scale - freq_spread_sum * freq_spread_sum;
    if (!quiet_run) printf("# Initial in cascade R: %lf Omega: %lf R_uni %lf\n", order_param, freq_spread, univ_order_param);
    adj_list[initial_edge_cut]->weight = 0;
    if (!quiet_run) printf("# Initial edge to cut: %d\n", initial_edge_cut);
    
    coupling_current = coupling_target;

    for (t = 1; t <= max_it_cascade; t++) {
        x_out = x + delta;

        compute_derivatives_with_failures(x, state_variables, derivatives, coupling_target, t);
        flag = r8_rkf45(compute_derivatives, num_nodes * 2, state_variables, derivatives, &x, x_out, &relerr, abserr, flag);
        // gbs_integrator(compute_derivatives, num_nodes * 2, state_variables, derivatives, &x, x_out, &relerr, abserr, &flag);
        x += delta;

        if(t == time_cascade[output_count-1]){
            z = 0;
            freq_spread_sum = 0;
            univ_order_param = 0;
            freq_spread_square = 0;
            total_coupling_weight = 0;

            if(do_community){
                memset(comm_sum_ord_param, 0, (num_communities + 1) * sizeof(double));
                memset(comm_sum_univ_ord_param, 0, (num_communities + 1) * sizeof(double));
                memset(comm_sum_freq_square_sum, 0, (num_communities + 1) * sizeof(double));
                memset(comm_sum_freq_sum, 0, (num_communities + 1) * sizeof(double));
                memset(comm_weights, 0, (num_communities + 1) * sizeof(double));

                for (i = 1; i <= num_nodes; i++) {
                    z += cexp(I * state_variables[i]);
                    node_frequency = state_variables[i + num_nodes];
                    freq_spread_square += node_frequency * node_frequency;
                    freq_spread_sum += node_frequency;

                    int comm_idx = nodes[i].community;
                    comm_sum_ord_param[comm_idx] += cexp(I * state_variables[i]);
                    comm_sum_freq_square_sum[comm_idx] += node_frequency * node_frequency;
                    comm_sum_freq_sum[comm_idx] += node_frequency;
                }

                for(i = 0; i < num_edges; i++){
                    double weight = adj_list[i]->weight;
                    int node1 = adj_list[i]->nodes[0];
                    int node2 = adj_list[i]->nodes[1];
                    univ_order_param += weight * cos(state_variables[node2] - state_variables[node1]);
                    total_coupling_weight += weight;
                
                    if (nodes[node1].community == nodes[node2].community) {
                        comm_sum_univ_ord_param[nodes[node1].community] += weight * cos(state_variables[node2] - state_variables[node1]);
                        comm_weights[nodes[node1].community] += weight;
                    }
                    
                }
                for (int i = 1; i <= num_communities; i++) {
                    if (comm_weights[i] > 0) {
                        comm_weights[i] = 1.0 / comm_weights[i];
                    }
                }
                
                for(i = 1; i <= num_communities; i++){
                    comm_ord_param_hist_cascade[output_count][i] = cabs(comm_sum_ord_param[i]) / (double)comm_count[i];
                    comm_univ_ord_param_hist_cascade[output_count][i] = comm_sum_univ_ord_param[i] * comm_weights[i];
                    comm_freq_spread_hist_cascade[output_count][i] = comm_sum_freq_square_sum[i] / (double)comm_count[i] - (comm_sum_freq_sum[i] / comm_count[i]) * (comm_sum_freq_sum[i] / comm_count[i]);
                }
                
            } else {

                for (i = 1; i <= num_nodes; i++) {
                    z += cexp(I * state_variables[i]);
                    node_frequency = state_variables[i + num_nodes];
                    freq_spread_square += node_frequency * node_frequency;
                    freq_spread_sum += node_frequency;
                }
                
                for(i = 0; i < num_edges; i++){
                    double weight = adj_list[i]->weight;
                    int node1 = adj_list[i]->nodes[0];
                    int node2 = adj_list[i]->nodes[1];
                    univ_order_param += weight * cos(state_variables[node2] - state_variables[node1]);
                    total_coupling_weight += weight;
                }
            }

            univ_order_param /= total_coupling_weight;

            order_param = cabs(z) * size_scale;
            freq_spread_sum = freq_spread_sum * size_scale;
            freq_spread = freq_spread_square * size_scale - freq_spread_sum * freq_spread_sum;

            ord_param_hist_cascade[output_count] = order_param;
            freq_spread_hist_cascade[output_count] = freq_spread;
            univ_ord_param_hist_cascade[output_count] = univ_order_param;
            failure_counts[output_count] = failure_count;

            if (!verbose) printf("%d %lf %lf %lf %lf\n",t, order_param, freq_spread, univ_order_param, coupling_current);
            
            output_count++;
            prev_order_param = order_param;
        }

    }
}

/***************************************************************************
 * @brief Allocates memory for system variables and tracking arrays.
 *
 * This function ensures that all necessary arrays and variables are
 * allocated dynamically before simulation begins.
 ***************************************************************************/
void pre_init_alloc() {
    /* ---- Allocate Memory for Logarithmic Time Binning ---- */
    time_cascade = (int *)malloc((max_time_cascade + 1) * sizeof(int));
    if (!time_cascade) {
        if (!time_cascade) {
            printf("Error: Memory allocation failed for time_cascade.\n");
            exit(EXIT_FAILURE);
        }
    }

    time_thermal = (int *)malloc((max_time_thermal + 1) * sizeof(int));
    if (!time_thermal) {
        if (!time_thermal) {
            printf("Error: Memory allocation failed for time_thermal.\n");
            exit(EXIT_FAILURE);
        }
    }

    lambda_thermal = (double *)malloc((max_it_thermal + 1) * sizeof(double));
    if (!lambda_thermal) {
        printf("Error: Memory allocation failed for lambda_thermal.\n");
        exit(EXIT_FAILURE);
    }

    ord_param_hist_th = (double*)malloc((max_time_thermal + 1) * sizeof(double));
    freq_spread_hist_th = (double*)malloc((max_time_thermal + 1) * sizeof(double));
    univ_ord_param_hist_th = (double*)malloc((max_time_thermal + 1) * sizeof(double));

    ord_param_hist_cascade = (double*)malloc((max_time_cascade + 1) * sizeof(double));
    freq_spread_hist_cascade = (double*)malloc((max_time_cascade + 1) * sizeof(double));
    univ_ord_param_hist_cascade = (double*)malloc((max_time_cascade + 1) * sizeof(double));
    failure_counts = (int*)malloc((max_time_cascade + 1) * sizeof(int));

    if (!ord_param_hist_th || !freq_spread_hist_th || !univ_ord_param_hist_th) {
        printf("Error: Memory allocation failed for history tracking.\n");
        exit(EXIT_FAILURE);
    }
    if (!ord_param_hist_cascade || !freq_spread_hist_cascade || !univ_ord_param_hist_cascade || !failure_counts) {
        printf("Error: Memory allocation failed for history tracking.\n");
        exit(EXIT_FAILURE);
    }
}

void post_init_alloc() {
    int i;

    last10Frequencies = (double**)malloc((num_nodes + 1) * sizeof(double*));
    for (i = 1; i <= num_nodes; i++) {
        last10Frequencies[i] = (double*)calloc((10 + 1), sizeof(double));
    }
    
    state_variables    = (double*)calloc((2 * num_nodes + 1), sizeof(double));
    derivatives        = (double*)calloc((2 * num_nodes + 1), sizeof(double));
    
    if (!state_variables || !derivatives) {
        printf("Error: Memory allocation failed for state vectors.\n");
        exit(EXIT_FAILURE);
    }

    local_ord_param = (double*)calloc((num_nodes + 1), sizeof(double));
    local_freq_spread = (double*)calloc((num_nodes + 1), sizeof(double));
    local_phase = (double*)calloc((num_nodes + 1), sizeof(double));

    if (do_community) {
        comm_ord_param_hist_therm = (double**)malloc((max_time_thermal + 1) * sizeof(double*));
        for (int i = 1; i <= max_time_thermal; i++) {
            comm_ord_param_hist_therm[i] = (double*)calloc(num_communities + 1, sizeof(double));
        }
        comm_freq_spread_hist_therm = (double**)malloc((max_time_thermal + 1) * sizeof(double*));
        for (int i = 1; i <= max_time_thermal; i++) {
            comm_freq_spread_hist_therm[i] = (double*)calloc(num_communities + 1, sizeof(double));
        }
        comm_univ_ord_param_hist_therm = (double**)malloc((max_time_thermal + 1) * sizeof(double*));
        for (int i = 1; i <= max_time_thermal; i++) {
            comm_univ_ord_param_hist_therm[i] = (double*)calloc(num_communities + 1, sizeof(double));
        }

        comm_ord_param_hist_cascade = (double**)malloc((max_time_cascade + 1) * sizeof(double*));
        for (int i = 1; i <= max_time_cascade; i++) {
            comm_ord_param_hist_cascade[i] = (double*)calloc(num_communities + 1, sizeof(double));
        }
        comm_freq_spread_hist_cascade = (double**)malloc((max_time_cascade + 1) * sizeof(double*));
        for (int i = 1; i <= max_time_cascade; i++) {
            comm_freq_spread_hist_cascade[i] = (double*)calloc(num_communities + 1, sizeof(double));
        }
        comm_univ_ord_param_hist_cascade = (double**)malloc((max_time_cascade + 1) * sizeof(double*));
        for (int i = 1; i <= max_time_cascade; i++) {
            comm_univ_ord_param_hist_cascade[i] = (double*)calloc(num_communities + 1, sizeof(double));
        }

        comm_sum_ord_param = (double*)calloc(num_communities + 1, sizeof(double));
        comm_sum_univ_ord_param = (double*)calloc(num_communities + 1, sizeof(double));
        comm_sum_freq_sum = (double*)calloc(num_communities + 1, sizeof(double));
        comm_sum_freq_square_sum = (double*)calloc(num_communities + 1, sizeof(double));
        comm_weights = (double*)calloc(num_communities + 1, sizeof(double));
    }

    failure_times = (double*)malloc((num_edges + 1) * sizeof(double));
    memset(failure_times, 0, (num_edges + 1) * sizeof(double));
    failure_nodes = (int**)malloc((num_edges + 1) * sizeof(int*));
    for (i = 1; i <= num_edges; i++) {
        failure_nodes[i] = (int*)malloc((2 + 1) * sizeof(int));
        memset(failure_nodes[i], 0, 3 * sizeof(int));
    }
}

void write_output() {
    int i, k;
    double power_flow;
    char filename[1000];
    double complex order_param_z;

    if (!time_thermal || !ord_param_hist_th || !freq_spread_hist_th || !univ_ord_param_hist_th) {
        printf("Error: Some thermal tracking arrays are uninitialized.\n");
        return;
    }
    if (!time_cascade || !ord_param_hist_cascade || !freq_spread_hist_cascade || !failure_counts || !univ_ord_param_hist_cascade) {
        printf("Error: Some cascade tracking arrays are uninitialized.\n");
        return;
    }

    sprintf(filename, "../results/%s/%s/a_%f/k_%f/simulation_results_l%f_a%f_t%f_sc%d_dp%f_id%ld.dat",
            io_dir, network_type, damping_coefficient, coupling_target, coupling_target, damping_coefficient, failure_threshold, scene, default_power, unique_id);
    
    printf("# Writing output to: %s\n", filename);
    FILE *fp = fopen(filename, "w");

    if (!fp) {
        printf("Error: Could not open output file %s\n", filename);
        exit(EXIT_FAILURE);
    }

    fprintf(fp, "### Kuramoto Simulation Results ###\n");
    fprintf(fp, "# Number of nodes: %d\n", num_nodes);
    fprintf(fp, "# Number of edges: %d\n", num_edges);
    fprintf(fp, "# Time increment: %lf\n", delta);
    fprintf(fp, "# Thermalization steps: %d\n", max_it_thermal);
    fprintf(fp, "# Cascade steps: %d\n", max_it_cascade);
    fprintf(fp, "# Initial Phase Randomness: %f\n", initial_phase_randomness);
    fprintf(fp, "# Initial line cut: %d\n", initial_edge_cut);
    fprintf(fp, "# Symmetric Network: %d\n", symmetric_network);
    fprintf(fp, "# Initial coupling: %lf\n", coupling_0);
    fprintf(fp, "# Target Coupling: %lf\n", coupling_target);
    fprintf(fp, "# Stages of coupling therm.: %d\n", num_stages);
    fprintf(fp, "# Damping coefficient: %lf\n", damping_coefficient);
    fprintf(fp, "# Failure threshold: %lf\n", failure_threshold);
    fprintf(fp, "\n");

    fprintf(fp, "### Order Parameter Evolution (Thermalization) ###\n");
    fprintf(fp, "# Time\tOrder Param\tFrequency Spread\tUniversal Order Param\n");

    for (i = 1; i <= max_time_thermal; i++) {
        fprintf(fp, "%d\t%lf\t%lf\t%lf",
            time_thermal[i-1],  
            ord_param_hist_th[i],
            freq_spread_hist_th[i],
            univ_ord_param_hist_th[i]); 

        if (do_community) {
            for (k = 1; k <= num_communities; k++) {
                fprintf(fp, "\t%lf\t%lf\t%lf", 
                        comm_ord_param_hist_therm[i][k],  
                        comm_freq_spread_hist_therm[i][k], 
                        comm_univ_ord_param_hist_therm[i][k]);
            }
        }

        fprintf(fp, "\n");
    }

    fprintf(fp, "\n### Order Parameter Evolution (Cascade) ###\n");
    fprintf(fp, "# Time\tOrder Param\tFrequency Spread\tUniversal Order Param\tFailures\n");

    for (i = 1; i <= max_time_cascade; i++) {
        fprintf(fp, "%d\t%lf\t%lf\t%lf\t%d",
            time_cascade[i-1],  
            ord_param_hist_cascade[i],
            freq_spread_hist_cascade[i],
            univ_ord_param_hist_cascade[i],
            failure_counts[i]);   

        if (do_community) {
            for (k = 1; k <= num_communities; k++) {
                fprintf(fp, "\t%lf\t%lf\t%lf", 
                        comm_ord_param_hist_cascade[i][k],  
                        comm_freq_spread_hist_cascade[i][k],
                        comm_univ_ord_param_hist_cascade[i][k]);
            }
        }

        fprintf(fp, "\n");
    }

    fprintf(fp, "\n### Failure Times & Nodes ###\n");
    fprintf(fp, "# Failure Time\tNode1\tNode2\n");
    for (i = 1; i <= num_edges; i++) {
        if (failure_times[i] == 0) continue;  
        fprintf(fp, "%lf\t%d\t%d\n",
            failure_times[i], failure_nodes[i][0], failure_nodes[i][1]);
    }
    fprintf(fp, "\n### Final Power Flow States ###\n");
    fprintf(fp, "# Node1\tNode2\tPower Flow\n");

    for (i = 0; i < num_edges; i++) {
        if (adj_list[i] == NULL) continue;  
        int node1 = adj_list[i]->nodes[0];
        int node2 = adj_list[i]->nodes[1];

        power_flow = coupling_target * adj_list[i]->weight * sin(state_variables[node2] - state_variables[node1]);
        fprintf(fp, "%d\t%d\t%lf\n", node1, node2, power_flow);
    }

    fprintf(fp, "\n### Final Local Kuramoto Phases & Frequencies Therm. ###\n");
    fprintf(fp, "# Node\tOrder Param\tFrequency Spread\tPhase\n");
    for (i = 1; i <= num_nodes; i++) {  
        fprintf(fp, "%d\t%lf\t%lf\t%lf\n", i, local_ord_param[i], local_freq_spread[i], local_phase[i]);
    }

    fprintf(fp, "\n### Final Local Kuramoto Phases & Frequencies Casc. ###\n");
    fprintf(fp, "# Node\tOrder Param\tFrequency Spread\tPhase\n");

    for (i = 1; i <= num_nodes; i++) {
        order_param_z = 0;
        double loc_freq_sum_square = 0, loc_freq_sum = 0;
        double loc_freq_spread = 0;

        for (k = 0; k < nodes[i].degree; k++) {
            int neighbor = nodes[i].neighbors[k];
            order_param_z += cexp(I * state_variables[neighbor]);
            loc_freq_sum += state_variables[neighbor + num_nodes];
            loc_freq_sum_square += state_variables[neighbor + num_nodes] * state_variables[neighbor + num_nodes];
        }

        double local_order_param = cabs(order_param_z) / (double)nodes[i].degree;
        loc_freq_sum /= (double)nodes[i].degree;
        loc_freq_spread = loc_freq_sum_square / (double)nodes[i].degree - loc_freq_sum * loc_freq_sum;

        fprintf(fp, "%d\t%lf\t%lf\t%lf\n", i, local_order_param, loc_freq_spread, state_variables[i]);
    }

    fprintf(fp, "\n### Failure Statistics ###\n");
    fprintf(fp, "# Total Failures: %d\n", failure_count);
    fprintf(fp, "# Final Universal Order Param: %lf\n", univ_ord_param_hist_cascade[max_it_cascade]);

    fprintf(fp, "\n### Last 10  frequencies ###\n");
    fprintf(fp, "\n#NodeID\t frequqncies\n");
    for (i = 1; i <= num_nodes; i++) {
        fprintf(fp, "%d", i);
        for (int j = 0; j < 10; j++) {
            fprintf(fp, "\t%lf", last10Frequencies[i][j]);
        }
        fprintf(fp, "\n");
    }
    
    fclose(fp);
    printf("# Output written successfully.\n");
}


