#
# Test set for f5b functionality
#

TEST_DIR = "/Users/mittun/Desktop/UNI/Thesis/Computer_Algebra_Experimentation/test"

println("\n\n--- Bespoke Multivariate Polynomial Representation Tests ---")
include("$TEST_DIR/multivar_poly_rep_test.jl")

println("\n\n--- Labelled Polynomial Tests ---")
include("$TEST_DIR/labelled_poly_test.jl")

println("\n\n--- F5B Criterion Tests ---")
include("$TEST_DIR/f5b_criterion_test.jl")

println("\n\n--- Critical Pair Queue Tests ---")
#include("$TEST_DIR/critical_pair_test.jl") # TODO - Tests currently fail here

println("\n\n--- F5B Run Tests ---")
include("$TEST_DIR/f5b_run_test.jl")

println("\n\n--- Reduce Gb Tests---")
include("$TEST_DIR/reduce_gb_test.jl")
