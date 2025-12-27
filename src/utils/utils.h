
#ifndef NUMERIC_COLUMNS_H
#define NUMERIC_COLUMNS_H

/*
 * ============================================================================
 *  Numeric Column Detection Utilities
 * ============================================================================
 *
 *  This header declares helper routines for detecting the number of numeric
 *  columns in text files containing purely numerical data.
 *
 *  The corresponding implementation assumes:
 *    - Input is a text file (not binary)
 *    - Rows are newline-separated
 *    - Columns are separated by whitespace, comma (,), semicolon (;),
 *      or pipe (|)
 *    - Fields are valid numeric literals accepted by strtod()
 *    - Comment lines begin with '#' (after optional leading whitespace)
 *    - The first non-comment, non-empty line contains all columns
 *
 * ============================================================================
 */

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Detect the number of numeric columns in a file stream.
 *
 * This function reads from the provided FILE stream until it encounters
 * the first non-empty, non-comment line. That line is assumed to contain
 * all numeric columns present in the file.
 *
 * Numeric fields may be integers, floating-point values, or scientific
 * notation and are parsed using strtod().
 *
 * @param fp  Open file stream positioned at the beginning of the file
 *
 * @return
 *   - > 0 : number of numeric columns detected
 *   - -1  : no data lines found in the file
 *   - -2  : malformed numeric data encountered
 *
 * @note
 *   The file position will be left immediately after the first data line.
 *   Callers should rewind(fp) if subsequent full parsing is required.
 */
int detect_numeric_columns(FILE *fp);

#ifdef __cplusplus
}
#endif

#endif /* NUMERIC_COLUMNS_H */

#ifndef SOLVER_FLAGS_H
#define SOLVER_FLAGS_H

/*
 * ============================================================================
 *  Solver Flag Handling Utilities
 * ============================================================================
 *
 *  This header declares a helper routine for handling and responding to
 *  status and error flags returned by numerical integration or solver
 *  routines.
 *
 *  The function interprets solver-specific flag values and applies
 *  corrective actions such as:
 *    - relaxing error tolerances
 *    - switching solver modes
 *    - retrying failed steps
 *    - terminating execution on fatal errors
 *
 *  This logic centralizes solver recovery behavior and avoids scattering
 *  flag-handling code throughout the integration loop.
 *
 * ============================================================================
 */

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Handle solver status flags and apply corrective actions.
 *
 * This function inspects the solver status flag pointed to by @p flag and
 * modifies solver parameters accordingly to allow the computation to
 * continue or terminate safely.
 *
 * The following flag values are interpreted:
 *
 *   - 3 : Relative error tolerance too small.
 *         RELERR is increased and the solver is instructed to continue.
 *
 *   - 5 : Solution vanished or underflow occurred.
 *         ABSERR is set to a small non-zero value and the solver retries
 *         in one-step mode.
 *
 *   - 6 : Requested accuracy could not be achieved.
 *         Both RELERR and ABSERR are relaxed.
 *
 *   - 7 : Integration routine is inefficient.
 *         Solver is instructed to continue, potentially switching modes.
 *
 *   - 8 : Invalid input parameters.
 *         Program execution is terminated.
 *
 * All other flag values are ignored.
 *
 * @param flag      Pointer to solver status flag (input/output)
 * @param relerr    Pointer to relative error tolerance (input/output)
 * @param abserr    Pointer to absolute error tolerance (input/output)
 * @param iteration Current solver iteration count (for diagnostics)
 *
 * @note
 *   - This function may modify the values pointed to by @p flag,
 *     @p relerr, and @p abserr.
 *   - On flag value 8, this function calls exit(1).
 *   - Diagnostic messages are printed to standard output.
 */
void handle_solver_flag(int *flag,
                        double *relerr,
                        double *abserr,
                        int iteration);

#ifdef __cplusplus
}
#endif

#endif /* SOLVER_FLAGS_H */
