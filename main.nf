#!/usr/bin/env nextflow

// This is a comment

/*
 * This is a block of comments
 */

// This is needed for activating the new DLS2
nextflow.enable.dsl=2


//Let's create a channel from string values

str = Channel.from('hello', 'hola', 'bonjour')

/*
* Let's print that channel using the operator view()
* https://www.nextflow.io/docs/latest/operator.html#view
*/
str.view()