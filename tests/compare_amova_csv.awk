#!/usr/bin/awk -f

BEGIN {
    FS = ","  
    diff_found = 0
}

NR == FNR {
    # read first file and store values
    key = $1","$2
    values[key] = sprintf("%.6f", $3)
    next
}

{
    # compare values with the second file
    key = $1","$2
    if (key in values) {
        val1 = values[key]
        val2 = sprintf("%.6f", $3)
        if (val1 != val2) {
            print "Difference found for " key ": " val1 " vs " val2""
            diff_found = 1
        }
    }
}

END {
    if (diff_found) {
        exit 1  
    }
}
