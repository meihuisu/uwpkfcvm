/**
 * @file test_api.c
 * @brief Bootstraps the test framework for the UWPKFCVM library.
 * @author - SCEC
 * @version 1.0
 *
 * Tests the UWPKFCVM library by loading it and executing the code as
 * UCVM would do it.
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "uwpkfcvm.h"

/**
 * Initializes and runs the test program. Tests link against the
 * static version of the library to prevent any dynamic loading
 * issues.
 *
 * @param argc The number of arguments.
 * @param argv The argument strings.
 * @return A zero value indicating success.
 */
int main(int argc, const char* argv[]) {

	// Declare the structures.
	uwpkfcvm_point_t pt;
	uwpkfcvm_properties_t ret;

	// Initialize the model.
        char *envstr=getenv("UCVM_INSTALL_PATH");
        if(envstr != NULL) {
            if (uwpkfcvm_init(envstr, "uwpkfcvm") != 0) {
                assert(1);
            }
            } else if (uwpkfcvm_init("..", "uwpkfcvm") != 0) {
                assert(1);
        }

	printf("Loaded the model successfully.\n");

	// Query a point.
        pt.longitude = -120.5050;
	pt.latitude = 35.960;
	pt.depth = 0;


	uwpkfcvm_query(&pt, &ret, 1);

	// vs/vp/rho are doubles ..
	assert(ret.vs > 0 && ret.vs == 1380.000000);
	assert(ret.vp > 0 && ret.vp == 2920.000000);
	assert(ret.rho > 0 && ret.rho == 2205.8792761710597);

	printf("Query was successful.\n");

	// Close the model.
	assert(uwpkfcvm_finalize() == 0);

	printf("Model closed successfully.\n");

	printf("\nALL UWPKFCVM TESTS PASSED\n");

	return 0;
}
