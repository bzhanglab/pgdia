nextflow.enable.dsl=2
include { sampleTissue; cohortTissue } from '../workflows/tissue_metadata'

workflow {
    assert sampleTissue([id: 'Cyno_Plex1_126_Adrenal__SRR123']) == 'Adrenal'
    assert sampleTissue([id: 'Cyno_Plex1_127N_Ovary__SRR123']) == 'Ovary'
    assert sampleTissue([id: 'Cyno_Plex2_127C_Prostate__ERR123']) == 'Prostate'
    assert sampleTissue([id: 'Cyno_Plex2_128N_GI_Ileum__DRR123']) == 'GI_Ileum'
    assert sampleTissue([id: 'Cyno_Plex12_131AB_Uterus_Cervix']) == 'Uterus_Cervix'
    assert sampleTissue([id: 'sample__RUN1', db_id: 'sample']) == 'sample'
    assert sampleTissue([id: 'PDAC_001']) == 'PDAC_001'
    assert sampleTissue([id: 'Cyno_Plex2_128N_GI_Ileum', tissue: 'Lymph node']) == 'Lymph_node'
    assert sampleTissue([id: 'Cyno_Plex2_128N_GI_Ileum', tissue: '  ']) == 'GI_Ileum'
    assert sampleTissue([id: 'Cyno_Plex2_128N_GI_Ileum', tissue: '']) == 'GI_Ileum'
    assert sampleTissue([id: 'Cyno_Plex2_128N_GI_Ileum', tissue: []]) == 'GI_Ileum'
    assert cohortTissue([[id: 'Cyno_Plex2_127N_Pancreas'], [id: 'Cyno_Plex2_130C_Heart']]) == 'pooled'
    assert cohortTissue([[id: 'run1', tissue: 'Pancreas'], [id: 'run2', tissue: 'Pancreas']]) == 'Pancreas'
    assert sampleTissue([id: 'all_rnaseq_runs', tissue: 'pooled']) == 'pooled'
    log.info 'Tissue metadata checks passed'
}
