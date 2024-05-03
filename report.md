### Introduction
Let matrix $A \in \mathbb{R}^{n\times n}$ is invertible, with the $spectral \ norm$ defined as $||A||_s = sup_{x\neq0}|Ax|/|x|$, the $condition\ number$ of $A$ is $c(A) = ||A||_s||A^{-1}||_s$. The $condition\ number$ is a measurement of the sensibility of the equation $Ax = b$ when right hand side is changed. If $c(A)$ is large, then $A$ is called $ill-conditioned$[[1](???)].

In [[2](???)], researcher restricted entries into set $\{0, 1\}$ or $\{-1, 1\}$, denoted by $\mathcal{A}_n^1$ and $\mathcal{A}_n^2$. With these conditions, many quantities are equivalent to the condition number. Let $A$ be a non-singular $(0, 1)$ matrix, $B = A^{-1} = (b_{ij})$, the following quantity is considered in [[2](???)]:
$$χ(A) = max_{i, j}|b_{ij}|.$$

In this report, we aim to construct a $(0, 1)$ matrix $C$ satisfied $χ(C) \geq2^{\frac{1}{2}n\log n-n(2+o(1))}$[[3](???)].
### 1 Generate Matrix $C$
#### 1.1 Generate Set $\Omega$
In order to generate matrix $A \in \mathcal{A}_n^2$, we need to generate set $\Omega$. Let $|\cdot|$ denotes cardinality and $\Delta$ denote symmetric different. Let $m\in \mathbb{Z}^{+}$, $n=2^m$, $\Omega = \{\alpha_0, \alpha_1, \alpha_2, ..., \alpha_n\}$ be a set of $n$ element such that $|\alpha_i| \leq |\alpha_{i+1}|$ and $|\alpha_i \Delta \alpha_{i+1}| \leq 2$. Let $\alpha_0=\{\varnothing\}$. We also have $\alpha_1=\{\varnothing\}$. Suppose we have $\Omega = \{\{\varnothing\}, \{\varnothing\}, \{1\}, \{2\}, ..., \{m\}\}$ initially. Considering for every time we only take out all the sets with maximum size in $\Omega$, and insert only one element inside by order. That is, for an existed set $\{1\}$, we insert $2, 3, ..., m$ by order.

With this way of insertion, we only need to consider if the conditions are met between the last old set and the first new set, and each of the continuous new sets comes from different old sets. For all the sets come from the same old set, the conditions are automatically met.
```
initially Omega = {{}, {}, {1}, {2}, ..., {m}}, we take out {1}, {2}, ..., {m}, 
and insert only one element. We need to make sure conditions are met between
{m} and {1, _}; {1, _} and {2, _}
```
For the first case above, the only element we can insert is $m$, which is the first element if we revert the ordered insertion. For the second case, let's say we have set $\alpha$ with size $k$ and $\beta$ with size $k-1$, if $\alpha \cup \beta = \alpha$, we can insert any element we wish; if $\alpha \cup \beta \neq \alpha$, we can only insert an element $r \in \alpha$. Since we always insert element in order or in reversed order, if the first element of the insertion list is not in $\alpha$, the first element of the reverted insertion list must be in $\alpha$ *(WHY???)*.

The following pseudo code shows the way generating $\Omega$: 
```
let Omega = {{}, {}, {1}, {2}, ..., {m}};
for i in {1, 2, ..., m - 1}: 
    let sets = sets in Omega with size equals i
    for set in sets:
        let avail_range = [max(set) + 1, ..., m]
        if avail_range[0] not in Omega[-1]:
            avail_range = avail_range.reverse()
        for element in avail_range:
            Omega.append({set..., element})  // set... means extend the set
```

when $m=4$, set $\Omega$ is shown below:
```
Omega = {set(), set(), {1}, {2}, {3}, {4}, 
         {1, 4}, {1, 3}, {1, 2}, {2, 4}, {2, 3}, 
         {3, 4}, {1, 3, 4}, {1, 2, 3}, {1, 2, 4}, 
         {2, 3, 4}, {1, 2, 3, 4}}
```

#### 1.2 Generate Matrix $A \in \mathcal{A}_n^2$
Shown in [[1](???)], with set $\Omega$ satisfying given conditions, we can generate matrix $A \in \mathcal{A}_n^2$ such that $χ(A)=2^{\frac{1}{2}n\log n-n(1+o(1))}$. Let $\mathcal{A}_n^2$ denotes the sets of invertible (-1, 1) matrices of order $n$. Let matrix $A \in \mathcal{A}_n^2$, with set $\Omega = \{\alpha_0, \alpha_1, \alpha_2, ..., \alpha_n\}$ satisfies $|\alpha_i| \leq |\alpha_{i+1}|$ and $|\alpha_i \Delta \alpha_{i+1}| \leq 2$, matrix $A$ can be constructed as follows[[1](???)].
For every $1 \leq i,j \leq n$:
$$a_{ij} = \begin{cases}
    -1,\;\alpha_j\bigcap(\alpha_{i-1}\bigcup\alpha_i)=\alpha_{i-1}\Delta\alpha_{i}\;and\;|\alpha_{i-1}\Delta\alpha_{i}|=2 \\
    (-1)^{|\alpha_{i-1}\bigcap\alpha_j| + 1},\;\alpha_{j}\bigcap(\alpha_{i-1}\bigcup\alpha_{i})\neq\varnothing\;but\;does\;not\;meet\;the\;condition\;above \\
    1,\;\alpha_j\bigcap(\alpha_{i-1}\bigcup\alpha_i)=\varnothing \\
\end{cases}.
$$

With the $\Omega$ shown in section 1.1 *(should be a link here)*, the matrix $A \in \mathcal{A}_n^2$ constructed is shown below:

$$
\left[\begin{smallmatrix}
    1 & 1 & 1 & 1 & 1 & 1 & 1 & 1 & 1 & 1 & 1 & 1 & 1 & 1 & 1 & 1 \\
    1 & -1 & 1  & 1  & 1  & -1 & -1 & -1 & 1  & 1  & 1  & -1 & -1 & -1 & 1  & -1 \\
    1 & 1  & -1 & 1  & 1  & 1  & 1  & -1 & -1 & -1 & 1  & 1  & -1 & -1 & -1 & -1 \\
    1 & 1  & 1  & -1 & 1  & 1  & -1 & 1  & 1  & -1 & -1 & -1 & -1 & 1  & -1 & -1 \\
    1 & 1  & 1  & 1  & -1 & -1 & 1  & 1  & -1 & 1  & -1 & -1 & 1  & -1 & -1 & -1 \\
    1 & -1 & 1  & 1  & 1  & 1  & -1 & -1 & 1  & 1  & 1  & 1  & -1 & 1  & 1  & 1  \\
    1 & 1  & 1  & -1 & 1  & -1 & 1  & 1  & 1  & -1 & -1 & -1 & 1  & -1 & -1 & -1 \\
    1 & 1  & -1 & 1  & 1  & 1  & -1 & 1  & -1 & -1 & 1  & -1 & -1 & 1  & -1 & -1 \\
    1 & 1  & 1  & 1  & -1 & -1 & 1  & -1 & 1  & 1  & -1 & -1 & -1 & -1 & 1  & -1 \\
    1 & 1  & 1  & -1 & 1  & 1  & -1 & 1  & -1 & 1  & -1 & -1 & 1  & -1 & -1 & -1 \\
    1 & 1  & 1  & 1  & -1 & -1 & 1  & 1  & -1 & -1 & 1  & 1  & -1 & -1 & -1 & -1 \\
    1 & -1 & 1  & 1  & 1  & 1  & 1  & -1 & 1  & 1  & -1 & -1 & 1  & 1  & -1 & -1 \\
    1 & 1  & -1 & 1  & 1  & -1 & -1 & 1  & -1 & 1  & -1 & 1  & -1 & -1 & -1 & 1  \\
    1 & 1  & 1  & 1  & -1 & 1  & -1 & -1 & 1  & -1 & -1 & -1 & 1  & -1 & -1 & 1  \\
    1 & 1  & 1  & -1 & 1  & -1 & -1 & -1 & -1 & 1  & 1  & -1 & -1 & 1  & -1 & 1  \\
    1 & -1 & 1  & 1  & 1  & 1  & 1  & 1  & -1 & -1 & -1 & -1 & -1 & -1 & 1  & 1  \\
\end{smallmatrix}\right]
$$

With the construction step above, it can be proved that $A=LQ$, where $Q$ is a $n$ by $n$ matrix given by $q_{ij} = (-1)^{|\alpha_i \bigcap \alpha_j|}$, and $Q$ is a symmetric Hadamard matrix, that is $Q^2=QQ^T=nI_n$; $L = AQ^{-1}$ is a lower triangular matrix.

With matrix $A$ shown above, matrix $Q$ and $L$ is shown below:
$$Q=
\left[\begin{smallmatrix}
    1 & 1 & 1 & 1 & 1 & 1 & 1 & 1 & 1 & 1 & 1 & 1 & 1 & 1 & 1 & 1 &  \\
    1 & -1 & 1 & 1 & 1 & -1 & -1 & -1 & 1 & 1 & 1 & -1 & -1 & -1 & 1 & -1 \\
    1 & 1 & -1 & 1 & 1 & 1 & 1 & -1 & -1 & -1 & 1 & 1 & -1 & -1 & -1 & -1 \\
    1 & 1 & 1 & -1 & 1 & 1 & -1 & 1 & 1 & -1 & -1 & -1 & -1 & 1 & -1 & -1 \\
    1 & 1 & 1 & 1 & -1 & -1 & 1 & 1 & -1 & 1 & -1 & -1 & 1 & -1 & -1 & -1 \\
    1 & -1 & 1 & 1 & -1 & 1 & -1 & -1 & -1 & 1 & -1 & 1 & -1 & 1 & -1 & 1 \\
    1 & -1 & 1 & -1 & 1 & -1 & 1 & -1 & 1 & -1 & -1 & 1 & 1 & -1 & -1 & 1 \\
    1 & -1 & -1 & 1 & 1 & -1 & -1 & 1 & -1 & -1 & 1 & -1 & 1 & 1 & -1 & 1 \\
    1 & 1 & -1 & 1 & -1 & -1 & 1 & -1 & 1 & -1 & -1 & -1 & -1 & 1 & 1 & 1 \\
    1 & 1 & -1 & -1 & 1 & 1 & -1 & -1 & -1 & 1 & -1 & -1 & 1 & -1 & 1 & 1 \\
    1 & 1 & 1 & -1 & -1 & -1 & -1 & 1 & -1 & -1 & 1 & 1 & -1 & -1 & 1 & 1 \\
    1 & -1 & 1 & -1 & -1 & 1 & 1 & -1 & -1 & -1 & 1 & -1 & 1 & 1 & 1 & -1 \\
    1 & -1 & -1 & -1 & 1 & -1 & 1 & 1 & -1 & 1 & -1 & 1 & -1 & 1 & 1 & -1 \\
    1 & -1 & -1 & 1 & -1 & 1 & -1 & 1 & 1 & -1 & -1 & 1 & 1 & -1 & 1 & -1 \\
    1 & 1 & -1 & -1 & -1 & -1 & -1 & -1 & 1 & 1 & 1 & 1 & 1 & 1 & -1 & -1 \\
    1 & -1 & -1 & -1 & -1 & 1 & 1 & 1 & 1 & 1 & 1 & -1 & -1 & -1 & -1 & 1 \\
\end{smallmatrix}\right]
$$

$$Q\times Q = 
\left[\begin{smallmatrix}
    16 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
    0 & 16 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
    0 & 0 & 16 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
    0 & 0 & 0 & 16&  0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
    0 & 0 & 0 & 0 & 16 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
    0 & 0 & 0 & 0 & 0 & 16 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
    0 & 0 & 0 & 0 & 0 & 0 & 16 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
    0 & 0 & 0 & 0 & 0 & 0 & 0 & 16 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
    0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 16 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
    0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 16 & 0 & 0 & 0 & 0 & 0 & 0 \\
    0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 16 & 0 & 0 & 0 & 0 & 0 \\
    0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 16&  0 & 0 & 0 & 0 \\
    0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 16 & 0 & 0 & 0 \\
    0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 16 & 0 & 0 \\
    0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 16 & 0 \\
    0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 16 \\
\end{smallmatrix}\right]
$$

$$L=
\left[\begin{smallmatrix}
    1.0 & 0.0 & 0.0 & 0.0 & 0.0 & 0.0 & 0.0 & 0.0 & 0.0 & 0.0 & 0.0 & 0.0 & 0.0 & 0.0 & 0.0 & 0.0 \\
    0.0 & 1.0 & 0.0 & 0.0 & 0.0 & 0.0 & 0.0 & 0.0 & 0.0 & 0.0 & 0.0 & 0.0 & 0.0 & 0.0 & 0.0 & 0.0 \\
    0.0 & 0.0 & 1.0 & 0.0 & 0.0 & 0.0 & 0.0 & 0.0 & 0.0 & 0.0 & 0.0 & 0.0 & 0.0 & 0.0 & 0.0 & 0.0 \\
    0.0 & 0.0 & 0.0 & 1.0 & 0.0 & 0.0 & 0.0 & 0.0 & 0.0 & 0.0 & 0.0 & 0.0 & 0.0 & 0.0 & 0.0 & 0.0 \\
    0.0 & 0.0 & 0.0 & 0.0 & 1.0 & 0.0 & 0.0 & 0.0 & 0.0 & 0.0 & 0.0 & 0.0 & 0.0 & 0.0 & 0.0 & 0.0 \\
    0.5 & 0.5 & 0.0 & 0.0 & -0.5 & 0.5 & 0.0 & 0.0 & 0.0 & 0.0 & 0.0 & 0.0 & 0.0 & 0.0 & 0.0 & 0.0 \\
    0.0 & 0.0 & 0.0 & 0.5 & 0.5 & -0.5 & 0.5 & 0.0 & 0.0 & 0.0 & 0.0 & 0.0 & 0.0 & 0.0 & 0.0 & 0.0 \\
    0.0 & 0.0 & 0.5 & 0.5 & 0.0 & 0.0 & -0.5 & 0.5 & 0.0 & 0.0 & 0.0 & 0.0 & 0.0 & 0.0 & 0.0 & 0.0 \\
    0.0 & 0.5 & 0.0 & 0.0 & 0.5 & 0.0 & 0.0 & -0.5 & 0.5 & 0.0 & 0.0 & 0.0 & 0.0 & 0.0 & 0.0 & 0.0 \\
    0.0 & 0.0 & 0.0 & 0.5 & 0.5 & 0.0 & 0.0 & 0.0 & -0.5 & 0.5 & 0.0 & 0.0 & 0.0 & 0.0 & 0.0 & 0.0 \\
    0.0 & 0.0 & 0.5 & 0.0 & 0.5 & 0.0 & 0.0 & 0.0 & 0.0 & -0.5 & 0.5 & 0.0 & 0.0 & 0.0 & 0.0 & 0.0 \\
    0.25 & 0.25 & 0.0 & 0.25 & 0.25 & 0.25 & 0.25 & 0.0 & 0.0 & 0.0 & -0.75 & 0.25 & 0.0 & 0.0 & 0.0 & 0.0 \\
    0.0 & 0.0 & 0.25 & 0.0 & 0.25 & 0.25 & 0.0 & 0.25 & 0.0 & 0.25 & 0.25 & -0.75 & 0.25 & 0.0 & 0.0 & 0.0 \\
    0.0 & 0.0 & 0.0 & 0.25 & 0.25 & 0.25 & 0.25 & 0.0 & 0.25 & 0.25 & 0.0 & 0.0 & -0.75 & 0.25 & 0.0 & 0.0 \\
    0.0 & 0.25 & 0.0 & 0.25 & 0.0 & 0.25 & 0.0 & 0.25 & 0.0 & 0.25 & 0.25 & 0.0 & 0.0 & -0.75 & 0.25 & 0.0 \\
    0.125 & 0.125 & 0.125 & 0.125 & 0.125 & 0.125 & 0.125 & 0.125 & 0.125 & 0.125 & 0.125 & 0.125 & 0.125 & 0.125 & -0.875 & 0.125\\
\end{smallmatrix}\right]
$$

#### 1.3 Generate Matrix $B \in \mathcal{A}_{n - 1}^1$
Let $\mathcal{A}_{n - 1}^1$ denotes the sets of invertible (0, 1) matrices of order $n$. Let matrix $B \in \mathcal{A}_{n - 1}^1$. Consider the map $\Phi$ which assigns to any matrix $B \in \mathcal{A}_{n - 1}^1$ a matrix $\Phi(B) \in \mathcal{A}_{n - 1}^1$ in the following way:
$$ \Phi(B) = \left(\begin{array}{cc} 
1 & 1_{n-1}\\
-1^T_{n-1} & 2B-J_{n-1}
\end{array}\right).
$$

Therefore, we have the following way to construct matrix $B \in \mathcal{A}_{n - 1}^1$ with $A=\{\alpha_{ij}\}\in \mathcal{A}_n^2$:

$$B = \frac{1}{2}(J_{n-1}+\{\alpha_{ij}\}_{2\le i\le n, 2\le j\le n}).$$

Notice that the $A \in \mathcal{A}_n^2$ we constructed above has it's first column as:
$$ A_2 = \Phi(B) = \left(\begin{array}{cc} 
1\\
1^T_{n-1}
\end{array}\right),
$$
we need a negative first column, which is said $\{\alpha_{ij}\}_{2\le i\le n, 2\le j\le n} \in \mathcal{A}_n^2$ multiplied with $-1$.  
Therefore we have the relation between matrix $B \in \mathcal{A}_{n - 1}^1$ and $A \in \mathcal{A}_n^2$ in the implementation shown as follows:    
$$B = \frac{1}{2}(J_{n-1}-\{\alpha_{ij}\}_{2\le i\le n, 2\le j\le n}).$$

With the matrix $A \in \mathcal{A}_n^2$ constructed above, a constructed matrix $B$ is shown below:

$$
\left[\begin{smallmatrix}
    1 & 0 & 0 & 0 & 1 & 1 & 1 & 0 & 0 & 0 & 1 & 1 & 1 & 0 & 1 \\
    0 & 1 & 0 & 0 & 0 & 0 & 1 & 1 & 1 & 0 & 0 & 1 & 1 & 1 & 1 \\
    0 & 0 & 1 & 0 & 0 & 1 & 0 & 0 & 1 & 1 & 1 & 1 & 0 & 1 & 1 \\
    0 & 0 & 0 & 1 & 1 & 0 & 0 & 1 & 0 & 1 & 1 & 0 & 1 & 1 & 1 \\
    1 & 0 & 0 & 0 & 0 & 1 & 1 & 0 & 0 & 0 & 0 & 1 & 0 & 0 & 0 \\
    0 & 0 & 1 & 0 & 1 & 0 & 0 & 0 & 1 & 1 & 1 & 0 & 1 & 1 & 1 \\
    0 & 1 & 0 & 0 & 0 & 1 & 0 & 1 & 1 & 0 & 1 & 1 & 0 & 1 & 1 \\
    0 & 0 & 0 & 1 & 1 & 0 & 1 & 0 & 0 & 1 & 1 & 1 & 1 & 0 & 1 \\
    0 & 0 & 1 & 0 & 0 & 1 & 0 & 1 & 0 & 1 & 1 & 0 & 1 & 1 & 1 \\
    0 & 0 & 0 & 1 & 1 & 0 & 0 & 1 & 1 & 0 & 0 & 1 & 1 & 1 & 1 \\
    1 & 0 & 0 & 0 & 0 & 0 & 1 & 0 & 0 & 1 & 1 & 0 & 0 & 1 & 1 \\
    0 & 1 & 0 & 0 & 1 & 1 & 0 & 1 & 0 & 1 & 0 & 1 & 1 & 1 & 0 \\
    0 & 0 & 0 & 1 & 0 & 1 & 1 & 0 & 1 & 1 & 1 & 0 & 1 & 1 & 0 \\
    0 & 0 & 1 & 0 & 1 & 1 & 1 & 1 & 0 & 0 & 1 & 1 & 0 & 1 & 0 \\
    1 & 0 & 0 & 0 & 0 & 0 & 0 & 1 & 1 & 1 & 1 & 1 & 1 & 0 & 0 \\
\end{smallmatrix}\right]
$$

#### 1.4 Generate and Verify Matrix $C$
Let $S$ and $T$ be two non-singular matrices of order $n_1$ and $n_2$. Define $S \diamond T$ ass follows:
$$
\left[\begin{smallmatrix}
    s_{11} & \dots & s_{1n_1} & 0 & \dots  & 0 \\
    s_{21} & \dots & s_{2n_1} & 0 & \dots  & 0 \\
    \vdots & \vdots & \vdots & \vdots & \vdots \\
    s_{n_11} & \dots & s_{n_1n_1} & 0 & \dots  & 0 \\
    0 & 0 \dots 0 & 1 & t_{11} & \dots & t_{1n_2} \\
    0 & 0 \dots 0 & 0 & t_{21} & \dots & t_{2n_2} \\
    \vdots & \vdots & \vdots & \vdots & \vdots \\
    0 & 0 \dots 0 & 0 & t_{n_21} & \dots & t_{n_2n_2} \\
\end{smallmatrix}\right]
$$

Consider the $(0, 1)$ matrix $C = A_1 \diamond (A_2 \diamond (. . . (A_{r−1} \diamond A_r))...)$. Let $M = C^{-1} = (m_{ij})$, $χ(C) = max_{i,j} |m_{ij}|$. Notice that $C$ is a spare $(0, 1)$ matrix, as shown in the article, $χ(C)$ has same order of magnitude as the condition number of matrix $C$, which can also be used for ill conditioned measurement. The following block shows some result of $χ(C)$ related to order $r$.

```
r = 2, order of C = 4, χ(C) = 1.0
r = 3, order of C = 11, χ(C) = 2.0
r = 4, order of C = 26, χ(C) = 260.0
r = 5, order of C = 57, χ(C) = 106491641548.6
```

We can see that as order $r$ grows, $χ(C)$ grows rapidly.

### 2 Complexity of Generating Matrix $C$
To generate set $\Omega$, everytime we only need to take out an generated set and insert one new element inside. Therefore, without cosidering the complexity of set insertion, the ideally time complexity is $O(n)$.

To generate every entry of matrix $A \in \mathcal{A}_n^2$, a visit of three continuous element in set $\Omega$ is necessary. Let $n$ be the order of matrix $A$, without considering the complexity of set accessing, the ideally time complexity is $O(n) + O(n^2) = O(n^2)$.

Generating matrix $B$ is a simple matrix operation, let $n - 1$ be the order of $B$, the time complexity is $O((n-1)^2) + O(n^2) = O(n^2)$.

In order to generate matrix $C$, we need a series of matrix $A_1$, $A_2$, ... , $A_r$. As shown above, for matrix $A_r$ with order $r$, the time complexity is $O(r^2)$. Therefore the time complexity generating these matrices is $\sum_{i=1}^{r}O(r^2) = O(n^3)$

Therefore, when $r$ grows, with an $O(n^3)$ time complexity, the time required to generate matrix $C$ grows rapidly.

> TODO: proofs on generation step? Code details? 