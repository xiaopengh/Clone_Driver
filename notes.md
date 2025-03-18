For a binomially distributed random variable $X \sim \text{Bin}(n, p_0)$, we use:

$$
P(X \leq k) = \sum_{i=0}^{k} \binom{n}{i} p_0^i (1 - p_0)^{n-i}
$$

For a one-sided test with "smaller" ($ H_1: p < p_0 $), the p-value is:

$$
P(X \leq k \mid H_0)
$$

For a one-sided test with "larger" ($ H_1: p > p_0 $), the p-value is:

$$
P(X \geq k \mid H_0) = 1 - P(X \leq k-1 \mid H_0)
$$

For a two-sided test ($ H_1: p \neq p_0 $), the p-value is:

$$
2 \times \min(P(X \leq k), P(X \geq k))
$$

