# METABAT
The project I did for my 2016 summer internship.

Use Evidence Accumulation Clustering on MetaBAT
Goal
∙
  Average out biases
∙
  Reduce variance
∙
  Avoid overfitting
∙
  Improve generalizability and robustness
Strategy
Apply one of Ensemble Clustering Method on MetaBAT --- Evidence Accumulation Clustering
∙
  Run N times MetaBAT on the same dataset.
∙
  Create a membership matrix which can record the binning membership for each contig
∙
  Record how many times each pair was clustered and create the co-association matrix
- For each run, group the row number together which share the same membership
- Label the pairs as 1 if they are in the same membership
- Get N sparse matrix and then average them co_assoc(i,j) = votes(i,j)/N
- Get the final similarity sparse matrix 
- The larger value shows contigs highly likely to be in the same cluster

In order to test if the algorithm can work well, run the algorithm on the subsample of the data first.
∙
  Subsample from the synthetic data which contains 195,601 contigs from 291 genomes
- Remove the genomes which have extreme number of contigs (>1000 or < 300)
- After removing, there are 250 genomes left in total
- Randomly choose 100 contigs from each genome ==> 25,000 contigs were finally chosen
- Sum the size of the contigs for each genome
- Get the summed size as range of 280k to 760k per genome
- Go through each genome and remove the contigs with larger size until the total size 
  of each genome can be between 280k and 330k. 
- Finally obtain 20466 contigs in total as the sample

Apply hierachical clustering on the similarity matrix of the sample data
∙
  Minimum spanning tree is one of the divisive hierarchical clustering methods, which start in one cluster and splits recursively as one moves down the hierachy
∙
  Use the number of times the two samples clustered together as the edge weights and obtain the weights from the co-association matrix
∙
  Apply MST and remove the edges with a weight smaller than the given threshold
Result
When using N = 200, which means running 200 times MetaBAT on the same subsample dataset, the following plot shows the F1 score for each-time running. We can see the best F1 score from the one-time run Metabat has the value about 0.85 and has 28 clusters in total.
