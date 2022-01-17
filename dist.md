# Analysis on K-mer frequency


## Size of K-mer

For DNA $k$-mer of size $k$, there should be in total $4^k$ different possible $k$-mers.

For a DNA sequence of length $L$, the total number of $k$-mers should be $L-k+1$. Assuming equal probability of each $k$-mer appearing, then for each particular $k$-mer, the expected number of times it appear in this sequence is 

$$\lambda=\frac{L-k+1}{4^k}$$

The number of times it appear in the sequence will follow a Poisson distribution with parameter $\lambda=(L-k+1)/4^k$.

## Reference