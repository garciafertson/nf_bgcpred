// main workflow assembly predict BGC using DeepBGC (and optional antismash)

nextflow.enable.dsl=2

include { BGCPRED } from './workflows/bgc_prediction'

//run assembly pipeline

workflow NF_BGC_PRED {
    BGCPRED ()
}


//     WORKFLOW: Execute a single named workflow for the pipeline

workflow {
    NF_BGC_PRED ()
}
