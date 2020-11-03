#!/usr/bin/env nextflow

// Some comment

/**
 * Some 
 * longer
 * comment
 */

// Creating a channel
numbers_ch = Channel.from(1,2,3)
strings_ch = Channel.from('a','b')

// Defining the process that is executed
process valuesToFile {

    input: 
    val nums from numbers_ch
    val strs from strings_ch
    
    output:
    file 'result.txt' into result_ch
    
    """
    echo $nums and $strs > result.txt
    """
}

// Using the operator view to inspect the result
result_ch.view{ "results file: $it"  }