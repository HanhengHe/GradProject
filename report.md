### 1 Generate Matrix $C$
#### 1.1 Generate Set $\Omega$
Let $|\cdot|$ denotes cardinality and $\Delta$ denote symmetric different. Let $m\in \mathbb{Z}^{+}$, $n=2^m$, $\Omega = \{\alpha_1, \alpha_2, ..., \alpha_n\}$ be a set of $n$ element such that $|\alpha_i| \leq |\alpha_{i+1}|$ and $|\alpha_i \Delta \alpha_{i+1}| \leq 2$. Let $\alpha_0=\{\varnothing\}$. The following pseudo code shows the way generating $\Omega$: 
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
Let $\mathcal{A}_n^2$ denotes the sets of invertible (-1, 1) matrices of order $n$. Let matrix $A \in \mathcal{A}_n^2$, with set $\Omega = \{\alpha_1, \alpha_2, ..., \alpha_n\}$ satisfies $|\alpha_i| \leq |\alpha_{i+1}|$ and $|\alpha_i \Delta \alpha_{i+1}| \leq 2$, matrix $A$ can be constructed as follows.
For every $1 \leq i,j \leq n$:
$$a_{ij} = \begin{cases}
    -1,\;\alpha_j\bigcap(\alpha_{i-1}\bigcup\alpha_i)=\alpha_{i-1}\Delta\alpha_{i}\;and\;|\alpha_{i-1}\Delta\alpha_{i}|=2 \\
    (-1)^{|\alpha_{i-1}\bigcap\alpha_j| + 1},\;\alpha_{j}\bigcap(\alpha_{i-1}\bigcup\alpha_{i})\neq\varnothing\;but\;does\;not\;meet\;the\;condition\;above \\
    1,\;\alpha_j\bigcap(\alpha_{i-1}\bigcup\alpha_i)=\varnothing \\
\end{cases}.
$$

With the $Omega$ shown in section 1.1, the matrix $A \in \mathcal{A}_n^2$ constructed is shown below:

```
┌─────────────────────────────────────────────────┐
 [1 1  1  1  1  1  1  1  1  1  1  1  1  1  1  1 ]
 [1 -1 1  1  1  -1 -1 -1 1  1  1  -1 -1 -1 1  -1]
 [1 1  -1 1  1  1  1  -1 -1 -1 1  1  -1 -1 -1 -1]
 [1 1  1  -1 1  1  -1 1  1  -1 -1 -1 -1 1  -1 -1]
 [1 1  1  1  -1 -1 1  1  -1 1  -1 -1 1  -1 -1 -1]
 [1 -1 1  1  1  1  -1 -1 1  1  1  1  -1 1  1  1 ]
 [1 1  1  -1 1  -1 1  1  1  -1 -1 -1 1  -1 -1 -1]
 [1 1  -1 1  1  1  -1 1  -1 -1 1  -1 -1 1  -1 -1]
 [1 1  1  1  -1 -1 1  -1 1  1  -1 -1 -1 -1 1  -1]
 [1 1  1  -1 1  1  -1 1  -1 1  -1 -1 1  -1 -1 -1]
 [1 1  1  1  -1 -1 1  1  -1 -1 1  1  -1 -1 -1 -1]
 [1 -1 1  1  1  1  1  -1 1  1  -1 -1 1  1  -1 -1]
 [1 1  -1 1  1  -1 -1 1  -1 1  -1 1  -1 -1 -1 1 ]
 [1 1  1  1  -1 1  -1 -1 1  -1 -1 -1 1  -1 -1 1 ]
 [1 1  1  -1 1  -1 -1 -1 -1 1  1  -1 -1 1  -1 1 ]
 [1 -1 1  1  1  1  1  1  -1 -1 -1 -1 -1 -1 1  1 ]
└─────────────────────────────────────────────────┘
```

#### 1.3 Generate Matrix $B \in \mathcal{A}_{n - 1}^1$
Let $\mathcal{A}_{n - 1}^1$ denotes the sets of invertible (0, 1) matrices of order $n$. Let matrix $B \in \mathcal{A}_{n - 1}^1$. Consider the map $\Phi$ which assigns to any matrix $B \in \mathcal{A}_{n - 1}^1$ a matrix $\Phi(B) \in \mathcal{A}_{n - 1}^1$ following way:
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
we need an adjustment, that is $\{\alpha_{ij}\}_{2\le i\le n, 2\le j\le n} \in \mathcal{A}_n^2$ multiplied with $-1$.  
Therefore we have the final relation between matrix $B \in \mathcal{A}_{n - 1}^1$ and $A \in \mathcal{A}_n^2$ in the implementation shown as follows:    
$$B = \frac{1}{2}(J_{n-1}-\{\alpha_{ij}\}_{2\le i\le n, 2\le j\le n}).$$

The matrix $B$ constructed is shown below:

```
┌───────────────────────────────┐
 [1 0 0 0 1 1 1 0 0 0 1 1 1 0 1]
 [0 1 0 0 0 0 1 1 1 0 0 1 1 1 1]
 [0 0 1 0 0 1 0 0 1 1 1 1 0 1 1]
 [0 0 0 1 1 0 0 1 0 1 1 0 1 1 1]
 [1 0 0 0 0 1 1 0 0 0 0 1 0 0 0]
 [0 0 1 0 1 0 0 0 1 1 1 0 1 1 1]
 [0 1 0 0 0 1 0 1 1 0 1 1 0 1 1]
 [0 0 0 1 1 0 1 0 0 1 1 1 1 0 1]
 [0 0 1 0 0 1 0 1 0 1 1 0 1 1 1]
 [0 0 0 1 1 0 0 1 1 0 0 1 1 1 1]
 [1 0 0 0 0 0 1 0 0 1 1 0 0 1 1]
 [0 1 0 0 1 1 0 1 0 1 0 1 1 1 0]
 [0 0 0 1 0 1 1 0 1 1 1 0 1 1 0]
 [0 0 1 0 1 1 1 1 0 0 1 1 0 1 0]
 [1 0 0 0 0 0 0 1 1 1 1 1 1 0 0]
└───────────────────────────────┘
```

#### Generate Matrix $C$
