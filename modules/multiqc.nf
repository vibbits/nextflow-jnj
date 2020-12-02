/*
*  multiqc module
*/
nextflow.enable.dsl=2

params.outdir = "$launchDir/results"

process multiqc {
    publishDir("$params.outdir/multiqc/", mode: 'copy', overwrite: true)
    label 'low'
    

    input:
    path (inputfiles)

    output:
    path "multiqc_report.html"					

    script:
    """
    multiqc .
    """
}




workflow MULTIQC {
    take: 
    input
    
    main:
		out =  qc(input)
    emit:
    	out
}
