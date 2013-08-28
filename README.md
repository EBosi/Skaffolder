SKAffolder is a generic scaffolding algorithm which approach relies on the presence of similar draft genomes used as reference.
The first module of SKA (Network construction) aligns similar sequences and evaluates synthenic regions to assess adjacency relationships between contigs, producing a scaffolding graph.
The second module (Network analysis), starting from the scaffolding graph, evaluates a local order/orientation of the contigs.

Skaffolder can accept as input a set of draft genomes, scaffolding graphs and orthologous sequences files.
The input files should be present in a Input folder (it can be set by user, by default is "Input", in the main directory), and the output files will be returned in a Output folder (it can be set by user, by default is "Output", in the main directory).
