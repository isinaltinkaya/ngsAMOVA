#include "dataStructs.h"
#include "em.h"
#include "metadata.h"










strArray* read_formula_str(const char* formula) {

    if (formula == NULL) {
        ERROR("No formula provided. Please provide a formula of the form y ~ x1/x2/.../xn via the --formula option.");
    }


    // pointer to the first character of formula string
    const char* p = formula;

    /// ------------------------------------------------------------
    /// get the first token - y (before the tilde ~)
    /// from the formula of form y ~ x1/x2/.../xn

    char* indtok = NULL;

    // skip until the tilde ~
    while (*p != '\0') {
        if (*p == '~') {
            ++p;
            // move the pointer to point to the character after the tilde after while loop
            break;
        }
        if (*p == '/') {
            // must not encounter any / before ~
            ERROR("Formula \"%s\" is not valid: Found '/' before '~'.", formula);
        }
        ++p;
    }

    // check if anything is left in the remaning formula string
    if (*p == '\0') {
        ERROR("Invalid formula \"%s\": No token found after '~'. Please provide a formula of the form y ~ x1/x2/.../xn.", formula);
    }

    char* token = strndup(formula, p - formula - 1);  // -1 to remove the tilde
    trimSpaces(token);

    indtok = strdup(token);

    strArray* formulaTokens = NULL;
    formulaTokens = strArray_init();

    IO::vprint(1, "Found the first token \"%s\" in formula \"%s\".\n", token, formula);

    /// ------------------------------------------------------------
    /// get remaining tokens - x(s) (after the tilde ~)
    /// if multiple, separated by /

    // point to the rest of the string
    const char* pstart = p;

    while (*p != '\0') {
        if (*p == '~') {
            fprintf(stderr, "\n[ERROR]\tFormula \"%s\" is not valid: Found more than one '~'. \n", formula);
            exit(1);
        }
        if (*p == '/') {
            FREE(token);
            token = strndup(pstart, p - pstart);
            trimSpaces(token);
            formulaTokens->add(token);
            IO::vprint(1, "Found new token \"%s\" in formula \"%s\".\n", token, formula);
            pstart = p + 1;
        }
        ++p;

        if (*p == '\0') {
            FREE(token);
            token = strndup(pstart, p - pstart);
            trimSpaces(token);
            formulaTokens->add(token);
            IO::vprint(1, "Found new token \"%s\" in formula \"%s\".\n", token, formula);
        }

    }

    formulaTokens->add(indtok);

    if (formulaTokens->len == 0) {
        ERROR("Formula \"%s\" is not valid: No tokens found.", formula);
    } else if (formulaTokens->len == 1) {
        ERROR("Formula \"%s\" is not valid: Found only one token (%s). Please provide a formula of the form y ~ x1/x2/.../xn.", formula, formula);
    }


    FREE(token);
    FREE(indtok);
    return (formulaTokens);
}

