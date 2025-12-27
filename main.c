/*****************************************************************
 *  Kuramoto Model Simulation                                   *
 *  This program simulates synchronization in a network using a *
 *  second-order Kuramoto model with community detection.       *
 *****************************************************************/

#include <omp.h>
#include <time.h>
#include <math.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <complex.h>
#include <sys/time.h>

#include "src/utils/utils.h"
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
int     seed;           
unsigned long long unique_id, state_id;        

/* ---- Network Parameters ---- */
int     num_nodes;        
int     num_edges;      
int     num_generators;  
int     initial_edge_cut; 
int     failure_count;    
int     num_cut_nodes;    
char    network_type[250];
int     random_self_freq;
int     force_freq_from_file;
int     load_state, save_state;
int     early_stopping_type;
int     early_stopping_time;
float   early_stopping_threshold;
int     saving_period;
int    shifted_node_id;
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
int     loadsave_it_th;   
int     loadsave_it_cs;   
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

Edge* adj_list;   

int symmetrize_network;

/* ---- Node Structure ---- */
typedef struct {
  int node_id;       
  int degree;        
  int community;     
  double power;
  double frequency;
  double gen_damping;
  double load_damping;
  double inertia;  
  int* neighbors;    
  int gen_index;
  int load_index;
} Node;

Node* nodes;

/* ---- Network Representation ---- */
double  *state_variables;
double  *derivatives;
double  *load_frequencies;
double  *coupling_sum;
double  *interaction;
double  *damping_array;
double  *inv_inertia_array;
/* ---- Output Files ---- */
char    in_dir[250];
char    out_dir[250];
char    slurm_job_id[250];
char    edges_filename[550];
char    nodes_filename[550];
char    state_filename[550];
char    state_fileids[250];

/* ---- Scene Identifier ---- */
int     scene;       

/* ---- Global History Tracking Arrays ---- */
double *ord_param_hist_th;
double *freq_spread_hist_th;
double *univ_ord_param_hist_th;

int max_save_loc;
double **local_ord_param;
double **local_freq_spread;
double **local_phase;
double **local_freq;

double **local_ord_param_out;
double **local_freq_spread_out;
double **local_phase_out;
double **local_freq_out;

double *ord_param_hist_cascade;
double *freq_spread_hist_cascade;
double *univ_ord_param_hist_cascade;

int *failure_counts;
double *failure_times;
int **failure_nodes;  
double **last10Frequencies; 
double **angle_diffs_hist;
double *last_order_params;

/* ---- Community Detection ---- */
int do_community; 
int num_communities;
int min_community_id;
int use_weight_runi;

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

/* ---- Normalizing DC ---- */
int norm_type; // 0: tanh, 1: tanh-like, 2: linear, 3: fcr
double norm_a;

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
    seed = 1976;
    num_edges = 0;
    num_nodes = 0;
    num_generators = 0;
    delta = 0.01;
    symmetrize_network = 1;
    random_self_freq = 1;
    force_freq_from_file = 0;
    precision_epsilon = 1.0e-12;
    max_it_thermal = 1000;
    max_it_cascade = 1000;
    load_state = 0;
    save_state = 0;
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
    num_communities = 0;
    min_community_id = 0;
    norm_type = 0;
    norm_a = 1.0;
    max_save_loc=10;
    use_weight_runi = 0;
    early_stopping_type = 1;
    early_stopping_time = 1000;
    early_stopping_threshold = 0.001;
    saving_period = 50000;
    shifted_node_id = 0;
    strcpy(in_dir, "hu387");
    strcpy(out_dir, "results");
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

        if (strcmp(argv[i], "-seed") == 0) seed = atoi(value);
        else if (strcmp(argv[i], "-init_rand") == 0) initial_phase_randomness = atof(value);
        else if (strcmp(argv[i], "-k") == 0) coupling_target = atof(value);
        else if (strcmp(argv[i], "-k0") == 0) coupling_0 = atof(value);
        else if (strcmp(argv[i], "-log_base") == 0) log_base = atof(value);
        else if (strcmp(argv[i], "-def_pow") == 0) default_power = atof(value);
        else if (strcmp(argv[i], "-it_th") == 0) max_it_thermal = atoi(value);
        else if (strcmp(argv[i], "-it_lc") == 0) max_it_cascade = atoi(value);
        else if (strcmp(argv[i], "-alpha") == 0) damping_coefficient = atof(value);
        else if (strcmp(argv[i], "-in_dir") == 0) strcpy(in_dir, value);
        else if (strcmp(argv[i], "-out_dir") == 0) strcpy(out_dir, value);
        else if (strcmp(argv[i], "-jid") == 0) strcpy(slurm_job_id, value);
        else if (strcmp(argv[i], "-delta") == 0) delta = atof(value);
        else if (strcmp(argv[i], "-stages") == 0) num_stages = atoi(value);
        else if (strcmp(argv[i], "-sym") == 0) symmetrize_network = atoi(value);
        else if (strcmp(argv[i], "-rsf") == 0) random_self_freq = atoi(value);
        else if (strcmp(argv[i], "-ff") == 0) force_freq_from_file = atoi(value);
        else if (strcmp(argv[i], "-ic") == 0) initial_edge_cut = atoi(value);
        else if (strcmp(argv[i], "-load_state") == 0) load_state = atoi(value);
        else if (strcmp(argv[i], "-save_state") == 0) save_state = atoi(value);
        else if (strcmp(argv[i], "-ls_it_cs") == 0) loadsave_it_cs = atoi(value);
        else if (strcmp(argv[i], "-ls_it_th") == 0) loadsave_it_th = atoi(value);
        else if (strcmp(argv[i], "-vb") == 0) verbose = atoi(value);
        else if (strcmp(argv[i], "-quiet") == 0) quiet_run = atoi(value);
        else if (strcmp(argv[i], "-thr") == 0) failure_threshold = atof(value);
        else if (strcmp(argv[i], "-scene") == 0) scene = atoi(value);
        else if (strcmp(argv[i], "-uwr") == 0) use_weight_runi = atoi(value);
        else if (strcmp(argv[i], "-id") == 0) unique_id = strtoull(value, NULL, 10);
        else if (strcmp(argv[i], "-nt") == 0) strcpy(network_type, value);
        else if (strcmp(argv[i], "-norm_typ") == 0) norm_type = atoi(value);
        else if (strcmp(argv[i], "-norm_a") == 0) norm_a = atof(value);
        else if (strcmp(argv[i], "-comms") == 0) do_community = atoi(value);
        else if (strcmp(argv[i], "-save_period") == 0) saving_period = atoi(value);
        else if (strcmp(argv[i], "-msl") == 0) max_save_loc = atoi(value); //TODO: use this at save last freqs
        else if (strcmp(argv[i], "-early") == 0) early_stopping_type = atoi(value);
        else if (strcmp(argv[i], "-earlyt") == 0) early_stopping_time = atoi(value);
        else if (strcmp(argv[i], "-earlyth") == 0) early_stopping_threshold = atof(value);
        else {
            printf("Error: Unknown argument %s\n", argv[i]);
            exit(1);
        }
    }
    if (max_it_thermal == 0) {
        max_it_thermal = 0;
    } else if (log_base <= 1.0) {
        max_time_thermal = max_it_thermal; 
    }else {
        int max_time_thermal_tmp = (int)(log((double)max_it_thermal)/log(log_base))+1;
        max_time_thermal = (max_time_thermal_tmp < max_it_thermal) ? max_time_thermal_tmp : max_it_thermal;
    }

    if (max_it_cascade == 0) {
        max_it_cascade = 0;
    } else if (log_base <= 1.0) {
        max_time_cascade = max_it_cascade; 
    } else {
        int max_time_cascade_tmp = (int)(log((double)max_it_cascade)/log(log_base))+1;
        max_time_cascade = (max_time_cascade_tmp < max_it_cascade) ? max_time_cascade_tmp : max_it_cascade;
    }
        
    i1 = gettimeofday(&tv, NULL);
    // if (!unique_id) unique_id = tv.tv_usec;
    
    unsigned long uunique_id = tv.tv_usec;
    // seed = uunique_id;
    // seed = -time(NULL);
    seed = (unsigned int)time(NULL) ^ getpid();

    // If running under SLURM, include job ID to make it unique
    char *slurm_job_id = getenv("SLURM_JOB_ID");
    if (slurm_job_id != NULL) {
        seed ^= (unsigned int)atoll(slurm_job_id);
    }

    if (load_state) {
        sprintf(state_fileids, "%s/%s/a_%f/k_%f/states/state_ids_l%f_a%f_t%f_sc%d_itth%d_itlc%d_initrand%f_normtyp%d_normint%f.dat",
            out_dir, network_type, damping_coefficient, coupling_target, coupling_target, damping_coefficient, failure_threshold, scene,
            loadsave_it_th, loadsave_it_cs, initial_phase_randomness, norm_type, norm_a);
        FILE *state_ids_file = fopen(state_fileids, "r");
        if (!state_ids_file) {
            fprintf(stderr, "Error: Could not open state IDs file %s\n", state_fileids);
            exit(EXIT_FAILURE);
        }
        printf("State file to read from: %s\n", state_fileids);
        unsigned long long id;
        size_t count = 0;
        while (fscanf(state_ids_file, "%llu", &id) == 1) {
            count++;
        }
    
        if (count == 0) {
            fprintf(stderr, "Error: No state IDs found in file\n");
            fclose(state_ids_file);
            exit(EXIT_FAILURE);
        }
    
        srand(seed);
        size_t random_index = rand() % count;
        printf("Randomly selecting state ID from %zu available IDs\n", count);
        printf("Random index: %zu\n", random_index);    
    
        rewind(state_ids_file);
        size_t current_index = 0;
        while (fscanf(state_ids_file, "%llu", &id) == 1) {
            if (current_index == random_index) {
                state_id = id;
                break;
            }
            current_index++;
        }
    
        fclose(state_ids_file);
    
        printf("Randomly selected state ID: %llu\n", state_id);
    }

    if (!verbose) printf("\n# Allocating memory pre init...\n");
    pre_init_alloc();
    
    if (!verbose) printf("# Initializing system...\n");
    initialize_system();

    if (!quiet_run) {
        printf("\n=== Simulation Parameters ===\n");
        printf("# System size (N)       : %d\n", num_nodes);
        printf("# Nr of generators (G)  : %d\n", num_generators);
        printf("# Nr of loads (L)       : %d\n", num_nodes - num_generators);
        printf("# Nr of edges (E)     : %d\n", num_edges);
        printf("# Thermalization steps  : %d\n", max_it_thermal);
        printf("# Cascade steps         : %d\n", max_it_cascade);
        printf("# Damping coefficient   : %lf\n", damping_coefficient);
        printf("# Target Coupling       : %lf\n", coupling_target);
        printf("# Initial Randomness    : %lf\n", initial_phase_randomness);
        printf("# Failure threshold     : %lf\n", failure_threshold);
        printf("# Cut edges             : %d\n", initial_edge_cut);
        printf("# Symmetrize network     : %d\n", symmetrize_network);
        printf("# Random self-freq      : %d\n", random_self_freq);
        printf("# Force freq from file  : %d\n", force_freq_from_file);
        printf("# Load sate             : %d\n", load_state);
        printf("# Save sate             : %d\n", save_state);
        printf("# Network type          : %s\n", network_type);
        printf("# Community analysis    : %d\n", do_community);
        if (do_community) {
            printf("# Nr of communities     : %d\n", num_communities);
        }
        printf("# Normalization type    : %d\n", norm_type);
        printf("# Linear Norm interval  : %lf\n", norm_a);
        printf("# Use weight runi       : %d\n", use_weight_runi);
        printf("# Logarithmic base      : %lf\n", log_base);
        printf("# Max save loc          : %d\n", max_save_loc);
        printf("# Early stopping type   : %d\n", early_stopping_type);
        printf("# Early stopping time   : %d\n", early_stopping_time);
        printf("# Early stopping threshold : %f\n", early_stopping_threshold);
        printf("# Saving period         : %d\n", saving_period);
        if (load_state) {
            printf("# Loaded state ID       : %llu\n", state_id);
        }
        printf("# ID                    : %lu\n", unique_id);
    }
    
    if (!verbose) printf("\n# Allocating memory post init...\n");
    post_init_alloc();

    if (!verbose) printf("\n# Generating network...\n");
    generate_network();

    if (!verbose) printf("\n# Starting simulation...\n");
    step();
    
    // if (!verbose) printf("\n# Writing output...\n");
    // write_output();
    
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
    printf("  -sym[bool]   : Symmetrize network edges\n");
    printf("  -rsf[bool]   : Random self-frequency assignment (0-centered Gaussian)\n");
    printf("  -ic[num]     : Initial edge ID to cut\n");
    printf("  -id[num]     : Unique ID\n");
    printf("  -load_state[bool] : Load state from previous run (0: no, 1: yes)\n");
    printf("  -save_state[bool] : Saves end state (0: no, 1: yes)\n");
    printf("  -early[num] : Early stopping type (1: order param, 2: angle diffs)\n");
    printf("  -earlyt[num] : Early stopping time (iterations after reaching target coupling)\n");
    printf("  -earlyth[num] : Early stopping threshold\n");
    printf("  -ls_it_th[num] : Iteration to load/save state from thermalization\n");
    printf("  -ls_it_cs[num] : Iteration to load/save state from cascade\n");
    printf("  -save_period[num] : Period for saving intermediate states\n");

    printf("\n---- Normalization DC ----\n");
    printf("  -norm_typ[num] : Normalization type (0: tanh, 1: tanh-like, 2: linear, 3: fcr)\n");
    printf("  -norm_a[num]   : Linear normalization interval\n");

    printf("\n---- Output Options ----\n");
    printf("  -in_dir[dir]  : input directory name\n");
    printf("  -out_dir[dir]  : output directory name\n");
    
    printf("\n---- Miscellaneous ----\n");
    printf("  -vb[num]     : Verbose output\n");
    printf("  -quiet[num]  : Suppress output\n");
    printf("  -uwr[bool]   : Use weight runi for coupling (default 0)\n");
    printf("  --h          : Display this help menu\n");

    printf("\nExample:\n");
    printf("  ./main -seed=1234 -k=0.5 -alpha=0.1 -it_th=1000 -it_lc=1000 -in_dir=input -out_dir=results -nt=base -scene=1\n\n");
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
    double base_frequency, weight, inertia, gen_damping, load_damping;
    int t, i, row, col, community_id, edge_type, force_freq_bool;
    force_freq_bool = 0;
    
    FILE *file;

    sprintf(edges_filename, "./networks/%s/%s/edges%d.data", in_dir, network_type, scene);

    file = fopen(edges_filename, "r");
    if (!file) {
        fprintf(stderr, "Error: Could not open edges file %s\n", edges_filename);
        exit(EXIT_FAILURE);
    }
    int min_node_id = 10;
    while (fscanf(file, "%d %d %lf %d\n", &row, &col, &weight, &edge_type) == 4) {
        num_edges++;
        if (row > num_nodes) num_nodes = row;
        if (col > num_nodes) num_nodes = col;
        if (row < min_node_id) min_node_id = row;
        if (col < min_node_id) min_node_id = col;
    }
    fclose(file);
    shifted_node_id = min_node_id;

    if (initial_edge_cut == -1) {
        srand(seed);
        initial_edge_cut = rand() % num_edges;
    }

    sprintf(nodes_filename, "./networks/%s/%s/nodes%d.data", in_dir, network_type, scene);
    file = fopen(nodes_filename, "r");

    if (!file) {
        printf("Error: Unable to open node data file: %s\n", nodes_filename);
        exit(EXIT_FAILURE);
    }

    int ncols = detect_numeric_columns(file);
    printf("# Detected %d columns in node data file.\n", ncols);
    if (ncols != 2 && ncols != 3 && ncols != 4 && ncols != 6 && ncols != 7) {
        fprintf(stderr,
            "Error: Unexpected number of columns (%d) in node data file: %s\n",
            ncols, nodes_filename);
            exit(EXIT_FAILURE);
        }
        
    rewind(file);

    const char *format = "%d %lf %lf %lf %lf %d %d";
    nodes = malloc(num_nodes * sizeof(Node));

    inv_inertia_array  = (double*)calloc(num_nodes, sizeof(double));
    damping_array      = (double*)calloc(num_nodes, sizeof(double));

    if (!nodes) {
        fprintf(stderr, "Error: Memory allocation failed for nodes.\n");
        exit(EXIT_FAILURE);
    }

    double total_sum = 0.0;
    double random_sum = 0.0;
    int random_count = 0;
    int gen_idx = 0;
    int load_idx = 0;

    for (i = 0; i < num_nodes; i++) {

        /* Defaults */
        int row = -1;
        double power = 0.0;
        double gen_damping = 0.0;
        double load_damping = 0.0;
        double inertia = 0.0;
        int community_id = -1;
        int force_freq_bool = 0;

        char line[8192];

        while (fgets(line, sizeof(line), file)) {
            char *p = line;

            while (isspace((unsigned char)*p)) p++;
            if (*p == '#' || *p == '\0' || *p == '\n')
                continue;

            int nread = sscanf(
                p,
                "%d %lf %d %d %lf %lf %lf",
                &row,
                &power,
                &community_id,
                &force_freq_bool,
                &load_damping,
                &gen_damping,
                &inertia
            );

            if ((ncols == 2 && nread < 2) ||
                (ncols == 3 && nread < 3) ||
                (ncols == 4 && nread < 4) ||
                (ncols == 6 && nread < 6) ||
                (ncols == 7 && nread < 7)) {
                printf("nread: %d\n", nread);
                fprintf(stderr,
                        "Error: Incorrect file format or insufficient data in %s\n",
                        nodes_filename);
                exit(EXIT_FAILURE);
            }

            break;
        }

        /* Remap fields depending on column count */
        if (ncols <= 2) {
            community_id     = 0;
            force_freq_bool  = 0;
            inertia          = 1.0;
            load_damping     = 1.0;
            gen_damping      = 0.0;
        }
        else if (ncols == 3) {
            force_freq_bool  = 0;
            inertia          = 1.0;
            load_damping     = 1.0;
            gen_damping      = 0.0;
        }
        else if (ncols == 4 || ncols == 5) {
            inertia          = 1.0;
            load_damping     = 1.0;
            gen_damping      = 0.0;
        }
        else if (ncols == 6) {
            gen_damping      = 0.0;
        }

        /* Bookkeeping */
        if (community_id > num_communities)
            num_communities = community_id;
        if (i==0)
            min_community_id = community_id;
        else if (community_id < min_community_id)
            min_community_id = community_id;
        row -= shifted_node_id;

        nodes[i].node_id      = row;
        nodes[i].degree       = 0;
        nodes[i].community    = community_id;
        nodes[i].neighbors    = NULL;
        nodes[i].gen_damping  = gen_damping;
        nodes[i].load_damping = load_damping;
        nodes[i].inertia      = inertia;
        nodes[i].power        = power;
        damping_array[i]      = damping_coefficient * (gen_damping + load_damping);
        if (inertia > 0) {
            nodes[i].gen_index = gen_idx++;
            nodes[i].load_index = -1;
            inv_inertia_array[i] = 1.0 / inertia;
            num_generators++;
        } else {
            nodes[i].gen_index = -1;
            nodes[i].load_index = load_idx++;
            inv_inertia_array[i] = 0.0;
        }

        /* Frequency assignment logic (unchanged) */
        nodes[i].frequency = 50.0 *(1.0 + generateRandomFloat(0, 1) * 0.05); // Base frequency with small Gaussian noise

        total_sum += nodes[i].power;
    }

    fclose(file);
    printf("# Min community ID: %d\n", min_community_id);
    for (i = 0; i < num_nodes; i++) {
        nodes[i].community -= min_community_id;
    }
    // num_communities = num_communities - min_community_id; // Adjust for zero-based indexing

    printf("\n# Node 0: Id: %d, Power: %lf, Frequency: %lf, Community: %d, Gen Damping: %lf, Load Damping: %lf, Inertia: %lf, Gen Index: %d",
        nodes[0].node_id, 
        nodes[0].power,
        nodes[0].frequency,
        nodes[0].community,
        nodes[0].gen_damping,
        nodes[0].load_damping,
        nodes[0].inertia,
        nodes[0].gen_index
    );
    printf("\n# Node 1000: Id: %d, Power: %lf, Frequency: %lf, Community: %d, Gen Damping: %lf, Load Damping: %lf, Inertia: %lf, Gen Index: %d",
        nodes[999].node_id, 
        nodes[999].power,
        nodes[999].frequency,
        nodes[999].community,
        nodes[999].gen_damping,
        nodes[999].load_damping,
        nodes[999].inertia,
        nodes[999].gen_index
    );

    if(do_community) {
        comm_count = (int*)calloc(num_communities, sizeof(int));
        for (i=0; i<num_nodes; i++) {
            comm_count[nodes[i].community]++;
        }
    }

    sample_logarithmically(max_it_thermal, max_time_thermal, log_base, 1, time_thermal);
    sample_logarithmically(max_it_cascade, max_time_cascade, log_base, 1, time_cascade);
    printf("# Time thermal: %d \n", time_thermal[0]);
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

    sprintf(edges_filename, "./networks/%s/%s/edges%d.data", in_dir, network_type, scene);
    file = fopen(edges_filename, "r");

    if (!file) {
        printf("Error: Unable to open edge file: %s\n", edges_filename);
        exit(EXIT_FAILURE);
    }

    adj_list = malloc(num_edges * sizeof(Edge));
    if (!adj_list) {
            fprintf(stderr,
                    "Error: Memory allocation failed for edge list\n");
            exit(EXIT_FAILURE);
        }

    if (!nodes || !adj_list) {
        printf("Error: Memory allocation failed for network storage.\n");
        exit(EXIT_FAILURE);
    }

    file = fopen(edges_filename, "r");

    if (!file) {
        printf("Error: Unable to reopen edge file: %s\n", edges_filename);
        exit(EXIT_FAILURE);
    }

    int ncols = detect_numeric_columns(file);

    if (ncols != 2 && ncols != 3 && ncols != 4) {
        fprintf(stderr,
                "Error: Unexpected number of columns (%d) in edge data file: %s\n",
                ncols, edges_filename);
        exit(EXIT_FAILURE);
    }

    rewind(file);

    total_coupling_weight = 0.0;
    int edge_index = 0;

    char line[8192];

    /* Superset format */
    const char *format = "%d %d %lf %d";

    while (fgets(line, sizeof(line), file)) {

        char *p = line;
        while (isspace((unsigned char)*p)) p++;

        if (*p == '#' || *p == '\0' || *p == '\n')
            continue;

        /* Defaults */
        int node1 = -1;
        int node2 = -1;
        double edge_weight = 1.0;
        int edge_type = 0;

        int nread = sscanf(p, format,
                        &node1,
                        &node2,
                        &edge_weight,
                        &edge_type);

        /* Validate */
        if ((ncols == 2 && nread < 2) ||
            (ncols == 3 && nread < 3) ||
            (ncols == 4 && nread < 4)) {

            fprintf(stderr,
                    "Error: Incorrect file format or insufficient data in %s\n",
                    edges_filename);
            exit(EXIT_FAILURE);
        }

        /* Remap missing columns */
        if (ncols == 2) {
            edge_weight = 1.0;
            edge_type   = 0;
        } else if (ncols == 3) {
            edge_type   = 0;
        }

        /* Remap node IDs */
        node1 -= shifted_node_id;
        node2 -= shifted_node_id;

        adj_list[edge_index].nodes[0]  = node1;
        adj_list[edge_index].nodes[1]  = node2;
        adj_list[edge_index].weight    = edge_weight;
        adj_list[edge_index].edge_type = edge_type;

        /* Update adjacency */
        nodes[node1].degree++;
        nodes[node2].degree++;

        nodes[node1].neighbors =
            realloc(nodes[node1].neighbors,
                    nodes[node1].degree * sizeof(int));
        nodes[node2].neighbors =
            realloc(nodes[node2].neighbors,
                    nodes[node2].degree * sizeof(int));

        if (!nodes[node1].neighbors || !nodes[node2].neighbors) {
            fprintf(stderr,
                    "Error: Memory reallocation failed for adjacency list\n");
            exit(EXIT_FAILURE);
        }

        nodes[node1].neighbors[nodes[node1].degree - 1] = node2;
        nodes[node2].neighbors[nodes[node2].degree - 1] = node1;

        /* Weight bookkeeping */
        if (use_weight_runi)
            total_coupling_weight += edge_weight;
        else
            total_coupling_weight += 1.0;

        edge_index++;
    }

    fclose(file);

    if (do_community) {
        memset(comm_weights, 0, (num_communities + 1) * sizeof(double));
        for (int i = 0; i < num_edges; i++) {
            int source = adj_list[i].nodes[0];
            int target = adj_list[i].nodes[1];

            if (nodes[source].community == nodes[target].community) {
                if (use_weight_runi) {
                    comm_weights[nodes[source].community] += adj_list[i].weight;
                } else {
                    comm_weights[nodes[source].community] += 1.0;
                }
            }
        }

        for (int i = 0; i < num_communities; i++) {
            if (comm_weights[i] > 0) {
                comm_weights[i] = 1.0 / comm_weights[i];
            }
        }
    }

    if (!quiet_run) printf("\n# Detected Nodes: %d, Edges: %d, Symmetrize %d, Total weight: %lf\n", num_nodes, num_edges, symmetrize_network, total_coupling_weight);
}

void record_line_failures(int time, double y[]) {
    int i;
    for (i = 0; i < num_edges; i++) {

        double weight = adj_list[i].weight;
        if (weight == 0.0) continue;

        if (adj_list[i].edge_type != 0) continue;

        int n1 = adj_list[i].nodes[0];
        int n2 = adj_list[i].nodes[1];

        double power_flow = sin(y[n2] - y[n1]);

        if (fabs(power_flow) > failure_threshold) {

            if (!verbose) {
                printf("# Line failure detected at time %d between nodes %d and %d\n",
                       time, n1, n2);
            }

            adj_list[i].weight = 0.0;

            failure_count++;
            failure_times[failure_count] = time;
            failure_nodes[failure_count][0] = n1;
            failure_nodes[failure_count][1] = n2;
        }
    }
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

void compute_derivatives(double x, double y[], double dydx[])
{
    int i;
    memset(coupling_sum, 0, num_nodes * sizeof(double));

    #pragma omp parallel for
    for (i = 0; i < num_edges; i++) {
        int n1 = adj_list[i].nodes[0];
        int n2 = adj_list[i].nodes[1];
        double w = adj_list[i].weight;
        if (w == 0.0) continue;

        double s = w * sin(y[n2] - y[n1]);
        interaction[i] = s;
        #pragma omp atomic
        coupling_sum[n1] += s;
        #pragma omp atomic
        coupling_sum[n2] -= s;
    }

    /* -------------------------------------------------
     * Node-local terms (no coupling yet)
     * ------------------------------------------------- */
    #pragma omp parallel for
    for (i = 0; i < num_nodes; i++) {
        if (nodes[i].inertia > 0.0) {
            /* generator */
            int k = nodes[i].gen_index;

            /* θ̇ = ω */
            dydx[i] = y[num_nodes + k];

            /* ω̇ = (p − Dω)/m */
            dydx[num_nodes + k] = (nodes[i].power - damping_array[i] * y[num_nodes + k]) * inv_inertia_array[i];

        } else {
            /* load: θ̇ = p / D */
            dydx[i] = nodes[i].power / damping_array[i];
            // load_frequencies[nodes[i].load_index] = dydx[i] + coupling_current * coupling_sum[i] / damping_array[i];
        }
    }
    /* -------------------------------------------------
     * Edge coupling
     * ------------------------------------------------- */   
    #pragma omp parallel for
    for (i = 0; i < num_edges; i++) {

        double weight = adj_list[i].weight;
        
        if (weight == 0.0) continue;
        
        int n1 = adj_list[i].nodes[0];
        int n2 = adj_list[i].nodes[1];

        int edge_type = adj_list[i].edge_type;

        /* --- frequencies for coupling --- */
        double omega1, omega2;

        if (nodes[n1].inertia > 0.0) {
            omega1 = y[num_nodes + nodes[n1].gen_index];
        } else {
            // omega1 = load_frequencies[nodes[n1].load_index];
            omega1 = 50.0;
        }

        if (nodes[n2].inertia > 0.0) {
            omega2 = y[num_nodes + nodes[n2].gen_index];
        } else {
            // omega2 = load_frequencies[nodes[n2].load_index];
            omega2 = 50.0;
        }

        double freq_diff = omega2 - omega1;
        double contrib1 = 0.0, contrib2 = 0.0;
        
        /* coupling law */
        if (edge_type == 1) {
            if (norm_type == 0) {
                contrib1 =  norm_tanh(freq_diff) * weight;
                contrib2 = -norm_tanh(freq_diff) * weight * (double)symmetrize_network;
            } else if (norm_type == 1) {
                contrib1 =  norm_softsign(freq_diff) * weight;
                contrib2 = -norm_softsign(freq_diff) * weight * (double)symmetrize_network;
            } else if (norm_type == 2) {
                contrib1 =  norm_linear(freq_diff, norm_a) * weight;
                contrib2 = -norm_linear(freq_diff, norm_a) * weight * (double)symmetrize_network;
            } else if (norm_type == 3) {
                contrib1 =  norm_fcr(freq_diff, 49.8, 50.2) * weight;
                contrib2 = -norm_fcr(freq_diff, 49.8, 50.2) * weight * (double)symmetrize_network;
            } else if (norm_type == 4) {
                contrib1 =  nodes[n1].power;
                contrib2 = -nodes[n2].power * (double)symmetrize_network;
            }  else if (norm_type == 5) {
                contrib1 = 0.0;
                contrib2 = 0.0;
            } else {
                printf("Error: Unknown normalization type %d\n", norm_type);
                exit(EXIT_FAILURE);
            }
        } else {
            contrib1 =  coupling_current * interaction[i];
            contrib2 = -coupling_current * interaction[i] * (double)symmetrize_network;
        }

        /* apply coupling */

        if (nodes[n1].inertia > 0.0) {
            int k1 = nodes[n1].gen_index;
            #pragma omp atomic
            dydx[num_nodes + k1] += contrib1 * inv_inertia_array[n1];
        } else {
            #pragma omp atomic
            dydx[n1] += contrib1 / damping_array[n1];
        }

        if (nodes[n2].inertia > 0.0) {
            int k2 = nodes[n2].gen_index;
            #pragma omp atomic
            dydx[num_nodes + k2] += contrib2 * inv_inertia_array[n2];
        } else {
            #pragma omp atomic
            dydx[n2] += contrib2 / damping_array[n2];
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
    int i, j, t;
    double complex z;
    int MAX_SAVE = 10;
    int stage_count = 0;
    int output_count = 0;
    int failure_count_previous = 0;
    double size_scale = 1.0 / (double) num_nodes;
    double x, x_out, order_param, prev_order_param, univ_order_param;
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

    // ================================
    // Initialization or state loading
    // ================================
    if (load_state) {
        sprintf(state_filename,
                "%s/%s/a_%f/k_%f/states/state_l%f_a%f_t%f_sc%d_itth%d_itlc%d_initrand%f_normtyp%d_normint%f_id%llu.dat",
                out_dir, network_type, damping_coefficient, coupling_target,
                coupling_target, damping_coefficient, failure_threshold, scene,
                loadsave_it_th, loadsave_it_cs,
                initial_phase_randomness, norm_type, norm_a, state_id);
        FILE *state_file = fopen(state_filename, "r");
        if (!state_file) {
            fprintf(stderr, "Error: Could not open state file %s\n", state_filename);
            exit(EXIT_FAILURE);
        }

        for (i = 0; i < num_nodes; i++) {
            if (fscanf(state_file, "%lf", &state_variables[i]) != 1) {
                fprintf(stderr, "Error: Failed to read state variable %d from file %s\n", i,
                        state_filename);
                fclose(state_file);
                exit(EXIT_FAILURE);
            }
        }
        for (i = 0; i < num_generators; i++) {
            if (fscanf(state_file, "%lf", &state_variables[i + num_nodes]) != 1) {
                fprintf(stderr, "Error: Failed to read frequency variable %d from file %s\n", i,
                        state_filename);
                fclose(state_file);
                exit(EXIT_FAILURE);
            }
        }
        fclose(state_file);
    } else {
        // #pragma omp parallel for reduction(+:z,freq_spread_sum,freq_spread_square)
        for (i = 0; i < num_nodes; i++) {
            double theta = initial_phase_randomness * PI2 * ran_lagged_fibonacci(&seed);
            double freq  = nodes[i].frequency;
            int k = nodes[i].gen_index;
            state_variables[i]           = theta;
            if (nodes[i].inertia > 0){
                state_variables[k + num_nodes] = freq;
                freq_spread_sum     += freq;
                freq_spread_square  += freq * freq;
            }

            z += cexp(I * theta);
        }
    }

    printf("First loaded state: %lf %lf\n", state_variables[1], state_variables[1 + num_nodes]);
    printf("First edge: %d %d %lf\n", adj_list[1].nodes[0], adj_list[1].nodes[1], adj_list[1].weight);

    #pragma omp parallel for reduction(+:univ_order_param)
    for (i = 0; i < num_edges; i++) {
        int node1 = adj_list[i].nodes[0];
        int node2 = adj_list[i].nodes[1];
        double weight = adj_list[i].weight;
        if (use_weight_runi)
            univ_order_param += weight * cos(state_variables[node2] - state_variables[node1]);
        else
            univ_order_param += cos(state_variables[node2] - state_variables[node1]);
    }

    univ_order_param /= total_coupling_weight;
    order_param       = cabs(z) * size_scale;
    prev_order_param  = 0;
    freq_spread_sum  /= (double)num_generators;
    freq_spread       = freq_spread_square / (double)num_generators - freq_spread_sum * freq_spread_sum;

    if (!quiet_run)
        printf("# Initial in therm. R: %lf Omega: %lf R_uni %lf\n", order_param, freq_spread,
               univ_order_param);

    // =============================
    // ===== THERMALIZATION ========
    // =============================

    if (!quiet_run) printf("\n# Stage %d, lambda %lf\n", stage_count + 1, lambda_thermal[1]);
    if (max_it_thermal != 0) {
        for (t = 0; t <= max_it_thermal; t++) {
            x_out = x + delta;
            coupling_current = lambda_thermal[t];

            if (t >= 1 && lambda_thermal[t] != lambda_thermal[t]) {
                stage_count++;
                if (!quiet_run)
                    printf("# Stage %d, lambda = %lf\n", stage_count + 1, coupling_current);
            }

            // Integrate one step
            compute_derivatives(x, state_variables, derivatives);
            flag = r8_rkf45(compute_derivatives, num_nodes * 2,
                            state_variables, derivatives,
                            &x, x_out, &relerr, abserr, flag);
            handle_solver_flag(&flag, &relerr, &abserr, t);
            x += delta;
            if ((early_stopping_type != 0 && t > early_stopping_time) || t == time_thermal[output_count] ) {
                // Reset accumulators
                z                = 0.0;
                freq_spread_sum  = 0.0;
                freq_spread_square = 0.0;
                univ_order_param = 0.0;
                // Community resets
                if (do_community) {
                    memset(comm_sum_ord_param, 0, (num_communities) * sizeof(double));
                    memset(comm_sum_univ_ord_param, 0, (num_communities) * sizeof(double));
                    memset(comm_sum_freq_square_sum, 0, (num_communities) * sizeof(double));
                    memset(comm_sum_freq_sum, 0, (num_communities) * sizeof(double));
                }

                // --- Node loop ---
                #pragma omp parallel for reduction(+:z,freq_spread_sum,freq_spread_square)
                for (i = 0; i < num_nodes; i++) {
                    double theta = state_variables[i];
                    double freq  = 0.0;
                    if (nodes[i].inertia > 0.0) {
                        int k = nodes[i].gen_index;
                        freq = state_variables[k + num_nodes];
                        last10Frequencies[k][output_count % MAX_SAVE] = freq;
                    }
                    // } else {
                    //     freq = load_frequencies[nodes[i].load_index];
                    // }

                    z += cexp(I * theta);
                    freq_spread_sum += freq;
                    freq_spread_square += freq * freq;
                    if (do_community) {
                        int comm = nodes[i].community;
                        #pragma omp atomic
                        comm_sum_ord_param[comm] += cexp(I * theta);
                        #pragma omp atomic
                        comm_sum_freq_sum[comm] += freq;
                        #pragma omp atomic
                        comm_sum_freq_square_sum[comm] += freq * freq;
                    }
                }

                order_param = cabs(z) * size_scale;
                last_order_params[t % 100] = order_param;

                
                // --- Edge loop ---
                #pragma omp parallel for reduction(+:univ_order_param)
                for (i = 0; i < num_edges; i++) {
                    int node1 = adj_list[i].nodes[0];
                    int node2 = adj_list[i].nodes[1];
                    double weight = adj_list[i].weight;

                    double dtheta = wrap_to_pi(state_variables[node2] - state_variables[node1]);
                    angle_diffs_hist[t % 100][i] = dtheta;

                    double cosval = cos(dtheta);
                    if (use_weight_runi) univ_order_param += weight * cosval;
                    else univ_order_param += cosval;

                    if (do_community && nodes[node1].community == nodes[node2].community) {
                        int comm = nodes[node1].community;
                        // printf("COMMUNITY %d EDGE BETWEEN %d AND %d\n", comm, node1, node2);
                        #pragma omp atomic
                        comm_sum_univ_ord_param[comm] += (use_weight_runi ? weight * cosval : cosval);
                    }
                }
            }

            // --- Output checkpoint ---
            if (t == time_thermal[output_count]) {
                univ_order_param /= total_coupling_weight;
                freq_spread_sum = freq_spread_sum / (double)num_generators;
                freq_spread = freq_spread_square / (double)num_generators - freq_spread_sum * freq_spread_sum;
                
                ord_param_hist_th[output_count]      = order_param;
                freq_spread_hist_th[output_count]    = freq_spread;
                univ_ord_param_hist_th[output_count] = univ_order_param;
                
                if (do_community) {
                    for (i = 0; i < num_communities; i++) {
                        comm_ord_param_hist_therm[output_count][i] =
                            cabs(comm_sum_ord_param[i]) / (double)comm_count[i];
                        comm_univ_ord_param_hist_therm[output_count][i] =
                            comm_sum_univ_ord_param[i] * comm_weights[i];
                        comm_freq_spread_hist_therm[output_count][i] =
                            comm_sum_freq_square_sum[i] / (double)comm_count[i] -
                            (comm_sum_freq_sum[i] / comm_count[i]) *
                            (comm_sum_freq_sum[i] / comm_count[i]);
                    }
                }
                
                // Local stats
                #pragma omp parallel for private(i)
                for (i = 0; i < num_nodes; i++) {
                    double complex order_param_z = 0.0;
                    double loc_freq_sum_square = 0.0, loc_freq_sum = 0.0;

                    for (int k = 0; k < nodes[i].degree; k++) {
                        int neighbor = nodes[i].neighbors[k];
                        order_param_z += cexp(I * state_variables[neighbor]);
                        double freq  = 0.0;
                        if (nodes[i].inertia > 0.0) {
                            int k = nodes[i].gen_index;
                            freq = state_variables[k + num_nodes];
                        }
                        
                        loc_freq_sum += freq;
                        loc_freq_sum_square += freq * freq;
                    }

                    int time_slot = output_count % (max_save_loc);
                    local_ord_param[i][time_slot] = cabs(order_param_z) / (double)nodes[i].degree;

                    double avg_freq = loc_freq_sum / (double)nodes[i].degree;
                    local_freq_spread[i][time_slot] =
                        loc_freq_sum_square / (double)nodes[i].degree - avg_freq * avg_freq;

                    local_phase[i][time_slot] = state_variables[i];
                    local_freq[i][time_slot]  = state_variables[num_nodes + i];
                }

                if (!verbose)
                    printf("%d %lf %lf %lf %lf\n", t, order_param, freq_spread, univ_order_param, coupling_current);

                if (t > early_stopping_time && coupling_current == coupling_target && early_stopping_type != 0) {
                    if (early_stopping_type == 1) {
                        double ord_param_diff = 0;
                        int count_nonzero = 0;
                        for (i = 0; i < 100; i++) {
                            if (last_order_params[i] == 0) break;
                            count_nonzero++;
                            ord_param_diff += fabs(last_order_params[i] - order_param);
                        }
                        ord_param_diff /= (count_nonzero > 0) ? count_nonzero : 1;

                        if (ord_param_diff < early_stopping_threshold) {
                            if (!quiet_run)
                                printf("# Early thermalization termination at step %d, avg order param change %lf\n",
                                       t, ord_param_diff);
                            break;
                        }
                    } else if (early_stopping_type == 2) {
                        double max_change = 0.0;
                        for (int j = 0; j < num_edges; j++) {
                            double minv = 1e9, maxv = -1e9;
                            for (int k = 0; k < 100; k++) {
                                double val = angle_diffs_hist[k][j];
                                minv = fmin(minv, val);
                                maxv = fmax(maxv, val);
                            }
                            if (maxv - minv > max_change) max_change = maxv - minv;
                        }

                        printf("Max change across history: %lf\n", max_change);

                        if (max_change < early_stopping_threshold) {
                            if (!quiet_run)
                                printf("# Early thermalization termination at step %d, max change %lf\n",
                                       t, max_change);
                            break;
                        }
                    }
                }

                output_count++;
                prev_order_param = order_param;
            }

            if ((t % saving_period == 0 && t != 0)|| t == max_it_thermal) {
                // ======================
                // Final averages
                // ======================
                #pragma omp parallel for private(i)
                for (i = 0; i < num_nodes; i++) {
                    double sum_ord = 0.0, sum_freq_spread = 0.0;
                    for (int j = 0; j < max_save_loc; j++) {
                        sum_ord         += local_ord_param[i][j];
                        sum_freq_spread += local_freq_spread[i][j];
                    }
                    double norm = (double)(max_save_loc);
                    local_ord_param_out[i][0]     = sum_ord / norm;
                    local_freq_spread_out[i][0]   = sum_freq_spread / norm;
                    local_phase_out[i][0]         = state_variables[i];
                    local_freq_out[i][0]          = state_variables[num_nodes + i];
                }
                if (!quiet_run) printf("# Writing output at step ... %d\n", t);
                write_output();
            }
        }
    }

    // ======================
    // Final averages
    // ======================
    #pragma omp parallel for private(i)
    for (i = 0; i < num_nodes; i++) {
        double sum_ord = 0.0, sum_freq_spread = 0.0;
        for (int j = 0; j < max_save_loc; j++) {
            sum_ord         += local_ord_param[i][j];
            sum_freq_spread += local_freq_spread[i][j];
        }
        double norm = (double)(max_save_loc);
        local_ord_param_out[i][0]     = sum_ord / norm;
        local_freq_spread_out[i][0]   = sum_freq_spread / norm;
        local_phase_out[i][0]         = state_variables[i];
        local_freq_out[i][0]          = state_variables[num_nodes + i];
    }
    
    /* ============================== */
    /* ===== CASCADE FAILURE ======== */
    /* ============================== */

    z = 0;
    x = 0.0;
    output_count = 0;
    failure_count = 0;
    freq_spread_sum = 0;
    univ_order_param = 0;
    freq_spread_square = 0;
    coupling_current = coupling_target;

    // clear local histories
    #pragma omp parallel for private(i)
    for (i = 0; i < num_nodes; i++) {
        memset(local_ord_param[i], 0, (max_save_loc) * sizeof(double));
        memset(local_freq_spread[i], 0, (max_save_loc) * sizeof(double));
        memset(local_phase[i], 0, (max_save_loc) * sizeof(double));
        memset(local_freq[i], 0, (max_save_loc) * sizeof(double));
    }

    if (!quiet_run) {
        printf("\n Casc time %d %d %d %d\n",
            time_cascade[0], time_cascade[1], time_cascade[10], time_cascade[100]);
    }

    // initial order parameter + spread
    #pragma omp parallel for reduction(+:z,freq_spread_sum,freq_spread_square)
    for (i = 0; i < num_nodes; i++) {
        double theta = state_variables[i];
        double freq  = 0.0;
        if (nodes[i].inertia > 0.0) {
            int k = nodes[i].gen_index;
            freq = state_variables[k + num_nodes];
        }

        z += cexp(I * theta);
        freq_spread_sum += freq;
        freq_spread_square += freq * freq;
    }

    #pragma omp parallel for reduction(+:univ_order_param)
    for (i = 0; i < num_edges; i++) {
        double weight = adj_list[i].weight;
        int n1 = adj_list[i].nodes[0];
        int n2 = adj_list[i].nodes[1];
        if (use_weight_runi)
            univ_order_param += weight * cos(state_variables[n2] - state_variables[n1]);
        else
            univ_order_param += cos(state_variables[n2] - state_variables[n1]);
    }

    prev_order_param = 0;
    order_param = cabs(z) * size_scale;
    univ_order_param /= total_coupling_weight;
    freq_spread_sum /= (double)num_generators;
    freq_spread = freq_spread_square / (double)num_generators - freq_spread_sum * freq_spread_sum;

    if (!quiet_run)
        printf("# Initial in cascade R: %lf Omega: %lf R_uni %lf\n",
            order_param, freq_spread, univ_order_param);

    adj_list[initial_edge_cut].weight = 0;
    if (!quiet_run) printf("# Initial edge to cut: %d\n", initial_edge_cut);

    if (max_it_cascade != 0) {
        for (t = 0; t <= max_it_cascade; t++) {
            x_out = x + delta;

            record_line_failures(t, state_variables);
            compute_derivatives(x, state_variables, derivatives);
            flag = r8_rkf45(compute_derivatives, num_nodes * 2,
                            state_variables, derivatives,
                            &x, x_out, &relerr, abserr, flag);
            handle_solver_flag(&flag, &relerr, &abserr, t);
            x += delta;

            if (t == time_cascade[output_count]) {
                // reset accumulators
                z = 0;
                freq_spread_sum = 0;
                freq_spread_square = 0;
                univ_order_param = 0;
                total_coupling_weight = 0;

                if (do_community) {
                    memset(comm_sum_ord_param, 0, (num_communities) * sizeof(double));
                    memset(comm_sum_univ_ord_param, 0, (num_communities) * sizeof(double));
                    memset(comm_sum_freq_square_sum, 0, (num_communities) * sizeof(double));
                    memset(comm_sum_freq_sum, 0, (num_communities) * sizeof(double));
                    memset(comm_weights, 0, (num_communities) * sizeof(double));
                }

                // --- Node loop ---
                #pragma omp parallel for reduction(+:z,freq_spread_sum,freq_spread_square)
                for (i = 0; i < num_nodes; i++) {
                    double theta = state_variables[i];
                    double freq  = 0.0;
                    if (nodes[i].inertia > 0.0) {
                        int k = nodes[i].gen_index;
                        freq = state_variables[k + num_nodes];
                    }

                    z += cexp(I * theta);
                    freq_spread_sum += freq;
                    freq_spread_square += freq * freq;

                    if (do_community) {
                        int comm = nodes[i].community;
                        #pragma omp atomic
                        comm_sum_ord_param[comm] += cexp(I * theta);
                        #pragma omp atomic
                        comm_sum_freq_sum[comm] += freq;
                        #pragma omp atomic
                        comm_sum_freq_square_sum[comm] += freq * freq;
                    }
                }

                // --- Edge loop ---
                #pragma omp parallel for reduction(+:univ_order_param,total_coupling_weight)
                for (i = 0; i < num_edges; i++) {
                    double weight = adj_list[i].weight;
                    int n1 = adj_list[i].nodes[0];
                    int n2 = adj_list[i].nodes[1];
                    double cosval = cos(state_variables[n2] - state_variables[n1]);

                    if (use_weight_runi) {
                        univ_order_param += weight * cosval;
                        total_coupling_weight += weight;
                    } else {
                        univ_order_param += cosval;
                        total_coupling_weight += 1.0;
                    }

                    if (do_community && nodes[n1].community == nodes[n2].community) {
                        int comm = nodes[n1].community;
                        if (use_weight_runi) {
                            #pragma omp atomic
                            comm_sum_univ_ord_param[comm] += weight * cosval;
                            #pragma omp atomic
                            comm_weights[comm] += weight;
                        } else {
                            #pragma omp atomic
                            comm_sum_univ_ord_param[comm] += cosval;
                            #pragma omp atomic
                            comm_weights[comm] += 1.0;
                        }
                    }
                }

                // finalize communities
                if (do_community) {
                    for (i = 0; i < num_communities; i++) {
                        if (comm_weights[i] > 0)
                            comm_weights[i] = 1.0 / comm_weights[i];

                        comm_ord_param_hist_cascade[output_count][i] =
                            cabs(comm_sum_ord_param[i]) / (double)comm_count[i];
                        comm_univ_ord_param_hist_cascade[output_count][i] =
                            comm_sum_univ_ord_param[i] * comm_weights[i];
                        comm_freq_spread_hist_cascade[output_count][i] =
                            comm_sum_freq_square_sum[i] / (double)comm_count[i] -
                            (comm_sum_freq_sum[i] / comm_count[i]) *
                            (comm_sum_freq_sum[i] / comm_count[i]);
                    }
                }

                // finalize global stats
                univ_order_param /= total_coupling_weight;
                order_param = cabs(z) * size_scale;
                freq_spread_sum /= (double)num_generators;
                freq_spread = freq_spread_square / (double)num_generators - freq_spread_sum * freq_spread_sum;

                ord_param_hist_cascade[output_count] = order_param;
                freq_spread_hist_cascade[output_count] = freq_spread;
                univ_ord_param_hist_cascade[output_count] = univ_order_param;
                failure_counts[output_count] = failure_count;

                // --- Local stats ---
                #pragma omp parallel for private(i)
                for (i = 0; i < num_nodes; i++) {
                    double complex order_param_z = 0;
                    double loc_freq_sum = 0.0, loc_freq_sum_sq = 0.0;

                    for (int k = 0; k < nodes[i].degree; k++) {
                        int neigh = nodes[i].neighbors[k];
                        order_param_z += cexp(I * state_variables[neigh]);
                        double freq  = 0.0;
                        if (nodes[i].inertia > 0.0) {
                            int k = nodes[i].gen_index;
                            freq = state_variables[k + num_nodes];
                        }
                        loc_freq_sum += freq;
                        loc_freq_sum_sq += freq * freq;
                    }

                    int time_slot = output_count % (max_save_loc);
                    local_ord_param[i][time_slot] = cabs(order_param_z) / (double)nodes[i].degree;

                    double avg_freq = loc_freq_sum / (double)nodes[i].degree;
                    local_freq_spread[i][time_slot] =
                        loc_freq_sum_sq / (double)nodes[i].degree - avg_freq * avg_freq;

                    local_phase[i][time_slot] = state_variables[i];
                    local_freq[i][time_slot] = state_variables[num_nodes + i];
                }

                if (!verbose)
                    printf("%d %lf %lf %lf %lf\n",
                        t, order_param, freq_spread, univ_order_param, coupling_current);

                output_count++;
                prev_order_param = order_param;
            }

                if ((t % saving_period == 0 && t != 0)|| t == max_it_thermal) {
                // ======================
                // Final averages
                // ======================
                #pragma omp parallel for private(i)
                for (i = 0; i < num_nodes; i++) {
                    double sum_ord = 0.0, sum_freq_spread = 0.0;
                    for (int j = 0; j < max_save_loc; j++) {
                        sum_ord += local_ord_param[i][j];
                        sum_freq_spread += local_freq_spread[i][j];
                    }
                    double norm = (double)max_save_loc;
                    local_ord_param_out[i][0] = sum_ord / norm;
                    local_freq_spread_out[i][0] = sum_freq_spread / norm;
                    local_phase_out[i][0] = state_variables[i];
                    local_freq_out[i][0] = state_variables[num_nodes + i];
                }

                if (!quiet_run) printf("# Writing output at step ... %d\n", t);
                write_output();
            }
        }
    }

    // final averages
    #pragma omp parallel for private(i)
    for (i = 0; i < num_nodes; i++) {
        double sum_ord = 0.0, sum_freq_spread = 0.0;
        for (int j = 0; j < max_save_loc; j++) {
            sum_ord += local_ord_param[i][j];
            sum_freq_spread += local_freq_spread[i][j];
        }
        double norm = (double)max_save_loc;
        local_ord_param_out[i][1] = sum_ord / norm;
        local_freq_spread_out[i][1] = sum_freq_spread / norm;
        local_phase_out[i][1] = state_variables[i];
        local_freq_out[i][1] = state_variables[num_nodes + i];
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
    time_cascade = (int *)malloc((max_time_cascade) * sizeof(int));
    if (!time_cascade) {
        if (!time_cascade) {
            printf("Error: Memory allocation failed for time_cascade.\n");
            exit(EXIT_FAILURE);
        }
    }

    time_thermal = (int *)malloc((max_time_thermal) * sizeof(int));
    if (!time_thermal) {
        if (!time_thermal) {
            printf("Error: Memory allocation failed for time_thermal.\n");
            exit(EXIT_FAILURE);
        }
    }
    
    lambda_thermal = (double *)malloc((max_it_thermal) * sizeof(double));
    if (!lambda_thermal) {
        printf("Error: Memory allocation failed for lambda_thermal.\n");
        exit(EXIT_FAILURE);
    }

    last_order_params = (double*)malloc((100) * sizeof(double));

    ord_param_hist_th = (double*)malloc((max_time_thermal) * sizeof(double));
    freq_spread_hist_th = (double*)malloc((max_time_thermal) * sizeof(double));
    univ_ord_param_hist_th = (double*)malloc((max_time_thermal) * sizeof(double));
    ord_param_hist_cascade = (double*)malloc((max_time_cascade) * sizeof(double));
    freq_spread_hist_cascade = (double*)malloc((max_time_cascade) * sizeof(double));
    univ_ord_param_hist_cascade = (double*)malloc((max_time_cascade) * sizeof(double));
    failure_counts = (int*)malloc((max_time_cascade) * sizeof(int));
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

    last10Frequencies = (double**)malloc((num_generators) * sizeof(double*));
    for (i = 0; i < num_generators; i++) {
        last10Frequencies[i] = (double*)calloc((10 + 1), sizeof(double));
    }
    
    state_variables    = (double*)calloc((num_nodes + num_generators), sizeof(double));
    derivatives        = (double*)calloc((num_nodes + num_generators), sizeof(double));
    load_frequencies   = (double*)calloc((num_nodes - num_generators), sizeof(double));
    coupling_sum       = (double*)calloc(num_nodes, sizeof(double));
    interaction        = (double*)calloc(num_edges, sizeof(double));

    
    if (!state_variables || !derivatives) {
        printf("Error: Memory allocation failed for state vectors.\n");
        exit(EXIT_FAILURE);
    }

    angle_diffs_hist = (double**)malloc((100) * sizeof(double*));
    for (i = 0; i < 100; i++) {
        angle_diffs_hist[i] = (double*)calloc((num_edges), sizeof(double));
    }

    local_ord_param = (double**)malloc((num_nodes) * sizeof(double*));
    for (i = 0; i < num_nodes; i++) {
        local_ord_param[i] = (double*)calloc((max_save_loc + 1), sizeof(double));
    }
    local_freq_spread = (double**)malloc((num_nodes) * sizeof(double*));
    for (i = 0; i < num_nodes; i++) {
        local_freq_spread[i] = (double*)calloc((max_save_loc + 1), sizeof(double));
    }
    local_phase = (double**)malloc((num_nodes) * sizeof(double*));
    for (i = 0; i < num_nodes; i++) {
        local_phase[i] = (double*)calloc((max_save_loc + 1), sizeof(double));
    }
    local_freq = (double**)malloc((num_nodes) * sizeof(double*));
    for (i = 0; i < num_nodes; i++) {
        local_freq[i] = (double*)calloc((max_save_loc + 1), sizeof(double));
    }

    local_ord_param_out = (double**)malloc((num_nodes) * sizeof(double*));
    for (i = 0; i < num_nodes; i++) {
        local_ord_param_out [i] = (double*)calloc((2), sizeof(double));
    }
    local_freq_spread_out  = (double**)malloc((num_nodes) * sizeof(double*));
    for (i = 0; i < num_nodes; i++) {
        local_freq_spread_out [i] = (double*)calloc((2), sizeof(double));
    }
    local_phase_out  = (double**)malloc((num_nodes) * sizeof(double*));
    for (i = 0; i < num_nodes; i++) {
        local_phase_out [i] = (double*)calloc((2), sizeof(double));
    }
    local_freq_out  = (double**)malloc((num_nodes) * sizeof(double*));
    for (i = 0; i < num_nodes; i++) {
        local_freq_out [i] = (double*)calloc((2), sizeof(double));
    }

    if (do_community) {
        comm_ord_param_hist_therm = (double**)malloc((max_time_thermal) * sizeof(double*));
        for (int i = 0; i < max_time_thermal; i++) {
            comm_ord_param_hist_therm[i] = (double*)calloc(num_communities, sizeof(double));
        }
        comm_freq_spread_hist_therm = (double**)malloc((max_time_thermal) * sizeof(double*));
        for (int i = 0; i < max_time_thermal; i++) {
            comm_freq_spread_hist_therm[i] = (double*)calloc(num_communities, sizeof(double));
        }
        comm_univ_ord_param_hist_therm = (double**)malloc((max_time_thermal) * sizeof(double*));
        for (int i = 0; i < max_time_thermal; i++) {
            comm_univ_ord_param_hist_therm[i] = (double*)calloc(num_communities, sizeof(double));
        }

        comm_ord_param_hist_cascade = (double**)malloc((max_time_cascade) * sizeof(double*));
        for (int i = 0; i < max_time_cascade; i++) {
            comm_ord_param_hist_cascade[i] = (double*)calloc(num_communities, sizeof(double));
        }
        comm_freq_spread_hist_cascade = (double**)malloc((max_time_cascade) * sizeof(double*));
        for (int i = 0; i < max_time_cascade; i++) {
            comm_freq_spread_hist_cascade[i] = (double*)calloc(num_communities, sizeof(double));
        }
        comm_univ_ord_param_hist_cascade = (double**)malloc((max_time_cascade) * sizeof(double*));
        for (int i = 0; i < max_time_cascade; i++) {
            comm_univ_ord_param_hist_cascade[i] = (double*)calloc(num_communities, sizeof(double));
        }

        comm_sum_ord_param = (double*)calloc(num_communities, sizeof(double));
        comm_sum_univ_ord_param = (double*)calloc(num_communities, sizeof(double));
        comm_sum_freq_sum = (double*)calloc(num_communities, sizeof(double));
        comm_sum_freq_square_sum = (double*)calloc(num_communities, sizeof(double));
        comm_weights = (double*)calloc(num_communities, sizeof(double));
    }

    failure_times = (double*)malloc((num_edges + 1) * sizeof(double));
    memset(failure_times, 0, (num_edges + 1) * sizeof(double));
    failure_nodes = (int**)malloc((num_edges + 1) * sizeof(int*));
    for (i = 1; i <= num_edges; i++) {
        failure_nodes[i] = (int*)malloc((3) * sizeof(int));
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

    sprintf(filename, "%s/%s/a_%f/k_%f/simulation_results_l%f_a%f_t%f_sc%d_itth%d_itlc%d_initrand%f_normtyp%d_normint%f_id%llu.dat",
            out_dir, network_type, damping_coefficient, coupling_target, coupling_target, damping_coefficient, failure_threshold, scene,
            max_it_thermal, max_it_cascade, initial_phase_randomness, norm_type, norm_a, unique_id);

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
    fprintf(fp, "# Symmetrize Network: %d\n", symmetrize_network);
    fprintf(fp, "# Initial coupling: %lf\n", coupling_0);
    fprintf(fp, "# Target Coupling: %lf\n", coupling_target);
    fprintf(fp, "# Stages of coupling therm.: %d\n", num_stages);
    fprintf(fp, "# Damping coefficient: %lf\n", damping_coefficient);
    fprintf(fp, "# Failure threshold: %lf\n", failure_threshold);
    fprintf(fp, "# Random self-freq      : %d\n", random_self_freq);
    fprintf(fp, "# Force freq from file  : %d\n", force_freq_from_file);
    fprintf(fp, "# Network type          : %s\n", network_type);
    fprintf(fp, "# Community analysis    : %d\n", do_community);
    fprintf(fp, "# Normalization type    : %d\n", norm_type);
    fprintf(fp, "# Linear Norm interval  : %d\n", norm_a);
    fprintf(fp, "# Use weight for R_uni  : %d\n", use_weight_runi);
    fprintf(fp, "# Unique ID             : %lu\n", unique_id);
    fprintf(fp, "# Maximum save loc.     : %d\n", max_save_loc);
    fprintf(fp, "# Initial edge cut      : %d\n", initial_edge_cut);
    fprintf(fp, "# Logarithmic base for time: %d\n", log_base);
    fprintf(fp, "# Save state            : %d\n", save_state);
    fprintf(fp, "# Load state            : %d\n", load_state);
    fprintf(fp, "# Early stopping type   : %d\n", early_stopping_type);
    fprintf(fp, "# Early stopping time   : %d\n", early_stopping_time);
    fprintf(fp, "# Early stopping threshold : %f\n", early_stopping_threshold);
    if (load_state) {
        fprintf(fp, "# Loaded state ID       : %llu\n", state_id);
        fprintf(fp, "# Load state itth       : %d\n", loadsave_it_th);
        fprintf(fp, "# Load state itlc       : %d\n", loadsave_it_cs);
    }
    fprintf(fp, "# Quiet run: %d\n", quiet_run);
    fprintf(fp, "# Verbose: %d\n", verbose);
    fprintf(fp, "\n");

    fprintf(fp, "### Order Parameter Evolution (Thermalization) ###\n");
    fprintf(fp, "# Time\tOrder Param\tFrequency Spread\tUniversal Order Param\n");

    for (i = 0; i < max_time_thermal; i++) {
        fprintf(fp, "%d\t%lf\t%lf\t%lf",
            time_thermal[i],  
            ord_param_hist_th[i],
            freq_spread_hist_th[i],
            univ_ord_param_hist_th[i]); 

        if (do_community) {
            for (k = 0; k < num_communities; k++) {
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

    for (i = 0; i < max_time_cascade; i++) {
        fprintf(fp, "%d\t%lf\t%lf\t%lf\t%d",
            time_cascade[i],  
            ord_param_hist_cascade[i],
            freq_spread_hist_cascade[i],
            univ_ord_param_hist_cascade[i],
            failure_counts[i]);   

        if (do_community) {
            for (k = 0; k < num_communities; k++) {
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
        // if (adj_list[i] == NULL) continue;  
        int node1 = adj_list[i].nodes[0];
        int node2 = adj_list[i].nodes[1];

        power_flow = coupling_target * adj_list[i].weight * sin(state_variables[node2] - state_variables[node1]);
        fprintf(fp, "%d\t%d\t%lf\n", node1, node2, power_flow);
    }

    fprintf(fp, "\n### Final Local Kuramoto Phases & Frequencies Therm. ###\n");
    fprintf(fp, "# Node\tOrder Param\tFrequency Spread\tPhase\tFrequency\n");
    for (i = 0; i < num_nodes; i++) {  
        fprintf(fp, "%d\t%lf\t%lf\t%lf\t%lf\n", i, local_ord_param_out[i][0], local_freq_spread_out[i][0], local_phase_out[i][0], local_freq_out[i][0]);
    }

    fprintf(fp, "\n### Final Local Kuramoto Phases & Frequencies Casc. ###\n");
    fprintf(fp, "# Node\tOrder Param\tFrequency Spread\tPhase\tFrequency\n");

    for (i = 0; i < num_nodes; i++) {  
        fprintf(fp, "%d\t%lf\t%lf\t%lf\t%lf\n", i, local_ord_param_out[i][1], local_freq_spread_out[i][1], local_phase_out[i][1], local_freq_out[i][1]);
    }

    fprintf(fp, "\n### Failure Statistics ###\n");
    fprintf(fp, "# Total Failures: %d\n", failure_count);
    fprintf(fp, "# Final Universal Order Param: %lf\n", univ_ord_param_hist_cascade[max_it_cascade]);

    fprintf(fp, "\n### Last 10  frequencies ###");
    fprintf(fp, "\n#NodeID\t frequencies\n");
    for (i = 0; i < num_generators; i++) {
        fprintf(fp, "%d", i);
        for (int j = 0; j < 10; j++) {
            fprintf(fp, "\t%lf", last10Frequencies[i][j]);
        }
        fprintf(fp, "\n");
    }
    
    fclose(fp);
    printf("# Output written successfully.\n");

    if (save_state) {
        sprintf(state_filename,  "%s/%s/a_%f/k_%f/states/state_l%f_a%f_t%f_sc%d_itth%d_itlc%d_initrand%f_normtyp%d_normint%f_id%llu.dat",
            out_dir, network_type, damping_coefficient, coupling_target, coupling_target, damping_coefficient, failure_threshold, scene, 
            max_it_thermal, max_it_cascade, initial_phase_randomness, norm_type, norm_a, unique_id);
        FILE *state_file = fopen(state_filename, "w");
        if (!state_file) {
            fprintf(stderr, "Error: Could not open state file %s for writing\n", state_filename);
            exit(EXIT_FAILURE);
        }
        
        for (i = 0; i < num_nodes; i++) {
            fprintf(state_file, "%lf\n", state_variables[i]);
        }
        for (i = 0; i < num_generators; i++) {
            fprintf(state_file, "%lf\n", state_variables[num_nodes + i]);
        }
        fclose(state_file);

        sprintf(state_fileids, "%s/%s/a_%f/k_%f/states/state_ids_l%f_a%f_t%f_sc%d_itth%d_itlc%d_initrand%f_normtyp%d_normint%f.dat",
            out_dir, network_type, damping_coefficient, coupling_target, coupling_target, damping_coefficient, failure_threshold, scene,
            max_it_thermal, max_it_cascade, initial_phase_randomness, norm_type, norm_a);
        FILE *id_file = fopen(state_fileids, "a");
        if (!id_file) {
            fprintf(stderr, "Error: Could not open id file %s for writing\n", state_fileids);
            exit(EXIT_FAILURE); 
        }
        fprintf(id_file, "%llu\n", unique_id);
        fclose(id_file);

        printf("# State saved successfully to %s and %s.\n", state_filename, state_fileids);
    }
}


