# Background Knowledge

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

1. **Suffix array interval (SA interval)**: If `W` is a substring of `X`, then the interval is \[*smallest index suffix with prefix `W`* in suffix array, index of *largest index suffix with prefix `W`*\]. These are denoted $[\underline{R}(W),\overline{R}(W)]$.

   Example: `X='GOOGOL`, `W=GO`, then the largest index suffix with prefix `W` is `GOOGOL` (`i=2,S[i]=0,B[i]=(googol)$`) while the smallest index is `GOL` (`i=1,S[i]=3,B[i]=(gol$go)o`), therefore the SA interval would be `[1,2]`.

1. **Backward Search**: For a character $a\in\Sigma$, and $W$ be a substring of $X$, we let

   $C(a):=$ the number of symbols in $X[0,n-2]$ that are lexicographically smaller than $a\in\Sigma$,

   $O(a,i):=$ the number of occurrences of $a$ in $B[0,i]$, then

   $$\underline{R}(aW)=C(a)+O(a,\underline{R}(W)-1)+1$$
   and
   $$\overline{R}(aW)=C(a)+O(a,\overline{R}(W))$$
   Then $aW$ is a substring of $X$ iff $\underline{R}(aW)\leq \overline{R}(aW)$. We can start from the end of $W$ and check one by one whether we have $\underline{R}(aW)\leq \overline{R}(aW)$ to test if $W$ is a substring.