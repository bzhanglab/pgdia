// Prefer explicit samplesheet metadata. Infer only the established Cyno naming
// convention; arbitrary sample identifiers are retained as fallback labels.
def sampleTissue(meta) {
    def tissue = meta.tissue ? meta.tissue.toString().trim() : ''
    if (!tissue) {
        tissue = (meta.db_id ?: meta.id ?: 'UNKNOWN').toString()
        if (!meta.db_id) {
            tissue = tissue.replaceFirst(/__(?:SRR|ERR|DRR)[0-9]+$/, '')
        }
        tissue = tissue.replaceFirst(/(?i)^Cyno_Plex[0-9]+_[0-9]+[A-Za-z]*_/, '')
    }
    // Keep Tissue= a single FASTA token, and quote-safe in process commands.
    return tissue.replaceAll(/[^A-Za-z0-9_.-]+/, '_') ?: 'UNKNOWN'
}

def cohortTissue(metas) {
    def tissues = metas.collect { sampleTissue(it) }.unique()
    // A merged assembly cannot establish a tissue of origin for each isoform.
    return tissues.size() == 1 ? tissues[0] : 'pooled'
}
