/*
 * ============================================================================
 *  Numeric Column Detection Utilities
 * ============================================================================
 *
 *  This module provides helper routines to determine the number of numeric
 *  columns in a text file containing purely numerical data.
 *
 *  Assumptions:
 *    - Input is a text file (not binary)
 *    - Rows are separated by newline characters
 *    - Columns are separated by standard delimiters:
 *        whitespace, comma (,), semicolon (;), or pipe (|)
 *    - Fields are valid numeric literals accepted by strtod()
 *      (integers, floating-point values, scientific notation)
 *    - Comment lines begin with '#' (after optional leading whitespace)
 *    - The first non-comment, non-empty line contains all columns present
 *      in the file
 *
 *  Under these assumptions, the column count is determined deterministically
 *  and without format-specific parsing logic.
 *
 *  Notes:
 *    - No memory allocation is performed
 *    - Parsing is locale-dependent via strtod()
 *    - This code is suitable for large files and streaming input
 *
 * ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

/**
 * @brief Count the number of numeric columns in a single line of text.
 *
 * This function scans a line and counts how many valid numeric values
 * it contains. Fields are assumed to be separated by one or more of the
 * following delimiters:
 *
 *   - whitespace (space, tab, etc.)
 *   - comma (,)
 *   - semicolon (;)
 *   - pipe (|)
 *
 * Numeric values are parsed using strtod(), allowing integers, floating-point
 * values, and scientific notation.
 *
 * @param line  Null-terminated string containing a single data line
 *
 * @return
 *   - >= 0 : number of numeric columns detected
 *   -  -1  : malformed numeric field encountered
 *
 * @note
 *   This function does not skip comment markers ('#') and expects to be
 *   called only on lines that contain data.
 */
static int count_numeric_columns_in_line(const char *line) {
    int cols = 0;
    const char *p = line;
    char *end;

    while (1) {
        /* Skip column delimiters */
        while (*p && (isspace((unsigned char)*p) ||
                      *p == ',' || *p == ';' || *p == '|'))
            p++;

        if (*p == '\0' || *p == '\n')
            break;

        /* Attempt to parse a numeric value */
        strtod(p, &end);
        if (end == p)
            return -1;  /* malformed numeric field */

        cols++;
        p = end;
    }

    return cols;
}

/**
 * @brief Detect the number of numeric columns in a file stream.
 *
 * This function reads lines from the provided FILE stream until it finds
 * the first non-empty, non-comment line. That line is assumed to contain
 * all columns present in the file.
 *
 * Comment lines are identified by a leading '#' character, optionally
 * preceded by whitespace.
 *
 * @param fp  Open file stream positioned at the beginning of the file
 *
 * @return
 *   - > 0 : number of numeric columns detected
 *   - -1  : no data lines found in the file
 *   - -2  : malformed numeric data encountered
 *
 * @note
 *   - The file position will be left after the first data line
 *   - The caller may rewind(fp) if subsequent full parsing is required
 */
int detect_numeric_columns(FILE *fp) {
    char line[8192];

    while (fgets(line, sizeof(line), fp)) {
        const char *p = line;

        /* Skip leading whitespace */
        while (isspace((unsigned char)*p))
            p++;

        /* Skip comments and empty lines */
        if (*p == '#' || *p == '\0' || *p == '\n')
            continue;

        int cols = count_numeric_columns_in_line(p);
        if (cols <= 0)
            return -2;

        /* First valid data line defines the column count */
        return cols;
    }

    /* No data lines found */
    return -1;
}

/*
* ============================================================================
*  Solver Flag Handling Utility
* ============================================================================
*
*  This module provides a helper routine to handle various solver return
*  flags by adjusting tolerances or modes as needed.
*
* ============================================================================
*/
void handle_solver_flag(int *flag, double *relerr, double *abserr, int iteration) {
    switch (*flag) {
        case 3:
            printf("# Warning: RELERR too small. Increasing RELERR. Iteration: %d\n", iteration);
            *flag = 2; // Reset flag to 2 to continue
            break;

        case 5:
            printf("# Warning: Solution vanished. Using non-zero ABSERR. Iteration: %d\n", iteration);
            if (*abserr == 0) *abserr = 1e-6; // Set ABSERR if zero
            *flag = 1; // One-step mode, try again
            break;

        case 6:
            printf("# Warning: Accuracy could not be achieved. Increasing tolerances. Iteration: %d\n", iteration);
            *relerr *= 2.0;  // Relax tolerance
            *abserr *= 2.0;  // Relax tolerance
            *flag = 2; // Reset flag to 2 to continue
            break;

        case 7:
            printf("# Warning: Integration routine inefficient. Trying one-step mode. Iteration: %d\n", iteration);
            *flag = 2; // Consider switching solver mode later if needed
            break;

        case 8:
            printf("# Error: Invalid input parameters. Exiting.\n");
            exit(1); // Exit program safely
            break;

        default:
            break; // No action needed for other flags
    }
}

