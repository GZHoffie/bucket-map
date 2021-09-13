# Notes

## Concepts

1. **Prefix trie**: The prefix trie of string `X` is a tree with edges being either `^` (start of the string) or a character in the string. Each path from a leaf to the root represent a substring of `X`.

1. **Burrows-Wheeler Transform**: Add a `$` sign at the end of the string, take all rotations of the string and sort them in lexicographic order. The last column of the matrix would be the BWT transform.

   Effects:
   - bringing together alike characters into runs.
   - Sort the characters according to the *right-context*, characters that come right after it. e.g. before the character `n` we tend to see vowels so vowels are sorted together.
   - BWM resembles the suffix array. We only sort the characters before the `$` sign (since when coming to the `$` sign, it is the smallest). So we are essentialy **sorting the suffix** of `X`.

     After sorting, we can just take the character before the sorted suffix to be the BWT.

     $$B[i]=\left\{\begin{array}{ll}X[S[i]-1]&S[i]>0\\\$&S[i]=0\end{array}\right.$$

   - **LF Mapping Character**: Relative order of same character remains the same for the last column and first column of BWM.

    Good reference: [cmu bioinfo lecture](https://www.cs.cmu.edu/~ckingsf/bioinfo-lectures/bwt.pdf)

1. **Suffix array interval (SA interval)**: If `W` is a substring of `X`, then the interval is \[*smallest index suffix with prefix `W`* in suffix array, index of *largest index suffix with prefix `W`*\]. These are denoted $[\underline{R}(W),\overline{R}(W)]$.

   Example: `X='GOOGOL`, `W=GO`, then the largest index suffix with prefix `W` is `GOOGOL` (`i=2,S[i]=0,B[i]=(googol)$`) while the smallest index is `GOL` (`i=1,S[i]=3,B[i]=(gol$go)o`), therefore the SA interval would be `[1,2]`.

1. **Backward Search**: For a character $a\in\Sigma$, and $W$ be a substring of $X$, we let

   $C(a):=$ the number of symbols in $X[0,n-2]$ that are lexicographically smaller than $a\in\Sigma$,

   $O(a,i):=$ the number of occurrences of $a$ in $B[0,i]$, then

   $$\underline{R}(aW)=C(a)+O(a,\underline{R}(W)-1)+1$$
   and
   $$\overline{R}(aW)=C(a)+O(a,\overline{R}(W))$$
   Then $aW$ is a substring of $X$ iff $\underline{R}(aW)\leq \overline{R}(aW)$. We can start from the end of $W$ and check one by one whether we have $\underline{R}(aW)\leq \overline{R}(aW)$ to test if $W$ is a substring.

## BWT Inexact Matching

**Idea**: Use a `D` array where `D[i]` stores the lower bound number of differences in `W[0,i]`.

**Algorithm**: use a recursive algorithm to calculate the SA interval of all similar matches. Kind of like a DFS on the prefix trie.

- `CalculateD(W)` function calculates the number of differences in `W[0, i]` and a substring of `X` (how many should be deleted to make it a substring of `X`). 

- `INEXRECUR(W, i, z, k, l)` recursively calculates the SA intervals of substrings that match `W[0, i]` with no more than `z` differences on the condition that suffix $W_{i+1}$ matches interval `[k, l]`.

**Features**:

1. Different penalties for mismatches and gaps.
2. heap-like data structure to store partial hits.
3. Iterative strategy: finds biggest interval by default.
4. Allow to set max differences.
5. save memory by calculating `O` and `S` arrays on the fly.
6. *What are SOLiD reads? (study later).*

## DNMFilter_Indel

**Goal**: to decrease false positive rate of de novo indels (rare, and high false positive rate).

**Method**:

- Training: Use state-of-the-art de novo indel detection methods; Use **synthetic**(simulated) and cross validated de novo indels as positive samples and random sampled indels as negative samples; refine alignment and extracts sequence features.

- Prediction: Use **gradient boosting classification model**, (so basically logical regression?)

## HAST

 Haplotype-resolved Assembly for Synthetic long reads using a Trio-binning strategy.

**Goal**: High haplotyping precision/recall

**Idea**: Use parental information to classify reads into maternal or paternal

**Tools**:

|Tools used|Purpose|
|--|--|
|SOAPfilter|filter low quality reads, duplicated reads and PE reads|
|TrioCanu (Canu)|classification of k-mers|
|Jellyfish|generate, count, output distinct k-mers in parental genomes|
|Merqury|k-mer analysis, haplotype evaluation|

**Methods**:

Hast steps:

1. Generate unshared k-mers between materal/pateral genomes.
2. *stLFR (What is that?)* sequencing of the child genome.
3. determine the origin of each long fragment in child's genome. Assign into 4 groups: paternal, ambiguous, maternal and homozygous. longer DNA fragments may be more easily categorized.
4. haplotype-resolved assemblies of both parental-inherited chromosomes.

Details:

1. Identifying maternal/paternal genome: have 1 copy of paternal k-mers and 2 copies of maternal k-mers and count the total frequencies.
2. limit parental k-mer library size based on the coverage distribution. Fit by a *mixture model of negative binomial model*. 
   k-mer with too large frequency (high-frequency duplications) and too small frequency (sequencing errors) are removed.
3. choose large `k` for small *collision rate* in large genomes.
4. To neutralize *inherent discrepancy* of parents, we
   - calculate the probability of the parent genome to be mapped by a DNA fragment.
   - Normailize to the parent k-mer library
   - Multiply by a correction factor of sex chromosome size variance.
5. binarize k-mer characters and use hash tables for k-mer database.


## Terms

|Vocabulary|Meaning|
|----------|-------|
|Phenotype|the observable physical properties of an organism|
|chromosomes|染色体|
|Diploid|二倍体 the number of each type of chromosome that an organism has|
|haplotype|单体型 a group of genes within an organism that was inherited together from a single parent|
|k-mer|substrings of length `k` contained within a biological sequence|
|heterozygosity|杂合性|
