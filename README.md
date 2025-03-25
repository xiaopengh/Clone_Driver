For a binomially distributed random variable $X \sim \text{Bin}(n, p_0)$, we use:

$$
P(X \leq k \mid p_0) = \sum_{i=0}^{k} \binom{n}{i} p_0^i (1 - p_0)^{n-i}
$$

For a one-sided test with alternative "smaller" i.e $H_1: p < p_0$, the p-value is:

$$
 P(X \leq k-1 \mid p_0) = 
 \text{inf} \{ P(X \leq k \mid p)  \mid p > p_0 \}
$$

For a one-sided test with alternative "larger" i.e $H_1: p > p_0$, the p-value is:

$$
P(X \geq k \mid p_0) = 1 - P(X \leq k-1 \mid p_0)
$$

For a two-sided test i.e $H_1: p \neq p_0$, the p-value is:

$$
2 \times \min(P(X \leq k), P(X \geq k))
$$

