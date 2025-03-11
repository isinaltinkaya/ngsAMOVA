/* ========================================================================== */
/// example general outline for header files:

/**
 * @file    example.h
 * @brief   header file for example.cpp
 * @details contains function declarations, macros, constants, and inline functions
 */

// #ifndef __EXAMPLE_H__
// #define __EXAMPLE_H__

/* SHARED ------------------------------------------------------------------- */
// #include "shared.h"  // shared definitions and global dependencies
/* END-OF-SHARED ---------------------------------------------------------------- */

/* INCLUDES ----------------------------------------------------------------- */
// #include <...>         // library headers
// #include "some_header.h"     // project-specific headers
/* END-OF-INCLUDES -------------------------------------------------------------- */

/* FORWARD-DECLARATIONS ----------------------------------------------------- */
// struct SomeStruct;    // forward declaration for a struct
/* END-OF-FORWARD-DECLARATIONS ------------------------------------------------ */

/* MACROS ------------------------------------------------------------------- */
// #define MAX_BUFFER_SIZE 1024
// #define ERROR_CODE -1
// #define SOME_MACRO_NAME(...) ...
/* END-OF-MACROS ---------------------------------------------------------------- */

/* CONSTANTS ---------------------------------------------------------------- */
// /**
//  * @brief    some_global_constant_lookup_table - lut for ...
//  * @details
//  *  maps ... to ...
//  */
// extern const int some_global_constant_lookup_table[256];  // global constant lookup table
/* END-OF-CONSTANTS ------------------------------------------------------------- */

/* ENUMS -------------------------------------------------------------------- */
// enum SomeEnum { ... };
/* END-OF-ENUMS ----------------------------------------------------------------- */

/* TYPEDEFS ----------------------------------------------------------------- */
// typedef unsigned int uint32_t;
/* END-OF-TYPEDEFS -------------------------------------------------------------- */

/* FUNCTIONS ---------------------------------------------------------------- */
// /**
//  * @brief square - returns the square of an integer
//  * @param x input integer
//  * @return squared value of x
//  */
// int square(int x);
/* END-OF-FUNCTIONS --------------------------------------------------------- */

/* TEMPLATES ---------------------------------------------------------------- */
// /**
//  * @brief returns the square of a value
//  * @tparam T type of the input
//  * @param value input value
//  * @return squared value of input
//  */
// template <typename T>
// T square(T value) {
//     return value * value;
// }
/* END-OF-TEMPLATES ------------------------------------------------------------- */

/* STATIC-INLINE ------------------------------------------------------------ */
// /**
//  * @brief adds two integers
//  * @param a first integer
//  * @param b second integer
//  * @return sum of a and b
//  */
// static inline int add(int a, int b) {
//     return a + b;
// }

// /**
//  * @brief subtracts two integers
//  * @param a first integer
//  * @param b second integer
//  * @return result of a - b
//  */
// static inline int subtract(int a, int b) {
//     return a - b;
// }
/* END-OF-STATIC-INLINE --------------------------------------------------------- */

// #endif  // __EXAMPLE_H__

/* ========================================================================== */

/// empty template to copy and paste for new header files

/**
 * @file    X.h
 * @brief   header file for X.cpp
 * @details contains function declarations, macros, constants, and inline functions
 */

// #ifndef __X_H__
// #define __X_H__

/* SHARED ------------------------------------------------------------------- */
// #include "shared.h"   // uncomment if needed
/* END-OF-SHARED ---------------------------------------------------------------- */

/* INCLUDES ----------------------------------------------------------------- */
/* END-OF-INCLUDES -------------------------------------------------------------- */

/* FORWARD-DECLARATIONS ----------------------------------------------------- */
/* END-OF-FORWARD-DECLARATIONS ------------------------------------------------ */

/* MACROS ------------------------------------------------------------------- */
/* END-OF-MACROS ---------------------------------------------------------------- */

/* CONSTANTS ---------------------------------------------------------------- */
/* END-OF-CONSTANTS ------------------------------------------------------------- */

/* ENUMS -------------------------------------------------------------------- */
/* END-OF-ENUMS ----------------------------------------------------------------- */

/* TYPEDEFS ----------------------------------------------------------------- */
/* END-OF-TYPEDEFS -------------------------------------------------------------- */

/* FUNCTIONS ---------------------------------------------------------------- */
/* END-OF-FUNCTIONS --------------------------------------------------------- */

/* TEMPLATES ---------------------------------------------------------------- */
/* END-OF-TEMPLATES ------------------------------------------------------------- */

/* STATIC-INLINE ------------------------------------------------------------ */
/* END-OF-STATIC-INLINE --------------------------------------------------------- */

// #endif  // __X_H__

/* ========================================================================== */