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

$$
\begin{bmatrix}
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
\end{bmatrix}
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
\begin{bmatrix}
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
\end{bmatrix}
$$

#### 1.4 Generate and Verify Matrix $C$
Let $S$ and $T$ be two non-singular matrices of order $n_1$ and $n_2$. Define $S \diamond T$ ass follows:
$$
\begin{bmatrix}
    s_{11} & \dots & s_{1n_1} & 0 & \dots  & 0 \\
    s_{21} & \dots & s_{2n_1} & 0 & \dots  & 0 \\
    \vdots & \vdots & \vdots & \vdots & \vdots \\
    s_{n_11} & \dots & s_{n_1n_1} & 0 & \dots  & 0 \\
    0 & 0 \dots 0 & 1 & t_{11} & \dots & t_{1n_2} \\
    0 & 0 \dots 0 & 0 & t_{21} & \dots & t_{2n_2} \\
    \vdots & \vdots & \vdots & \vdots & \vdots \\
    0 & 0 \dots 0 & 0 & t_{n_21} & \dots & t_{n_2n_2} \\
\end{bmatrix}
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