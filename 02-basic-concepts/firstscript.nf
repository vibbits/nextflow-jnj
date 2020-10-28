#!/usr/bin/env nextflow

// Some comment

/**
 * Some 
 * longer
 * comment
 */

numbers_ch = Channel.from(1,2,3)
strings_ch = Channel.from('a','b')

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

result_ch.view{ "results file: $it"  }