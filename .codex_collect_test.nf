nextflow.enable.dsl=2

workflow {
  Channel.of([1,2],[3,4])
    .collect()
    .view { "collect([1,2],[3,4]) -> ${it} (class=${it.getClass().name})" }

  Channel.of([1],[2],[3])
    .collect()
    .view { "collect([1],[2],[3]) -> ${it}" }

  Channel.of(1,2,3)
    .collect()
    .view { "collect(1,2,3) -> ${it}" }

  Channel.of([1,2],[3,4])
    .map { [it] }
    .collect()
    .view { "map{[it]}.collect() on list items -> ${it}" }
}
