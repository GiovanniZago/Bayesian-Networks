<!-- # R - K2 algorithm -->

<h1 align="center">Advanced Statistics - PoD<br>AY 2022/2023 <br> University of Padua</h1>

<p align="center">
  <img src="https://user-images.githubusercontent.com/62724611/166108149-7629a341-bbca-4a3e-8195-67f469a0cc08.png" height="150"/>
   
  <img src="https://user-images.githubusercontent.com/62724611/166108076-98afe0b7-802c-4970-a2d5-bbb997da759c.png" height="150"/>
</p>

<h3 align="center"><b>Group members:</b> E. Sarte, R. Tancredi, G. Zago<br></h3>

# Bayesian Networks

A **Bayesian Network**, known also as *belief network*, is a probabilistic graphical model that represents the relationships between given variables as a *directed acyclic graph* - **DAG**s **[1]**.

Each **node** in the graph corresponds to one of the variables $A, B, C, D \dots$ in the survey and the *direct* and *indirect* relations are expressed via **arcs** between pairs of variables: $A\to E$ means that $E$, the *child* node, depends on $A$, the *parent* **[2]**.

To create and manipulating BN DAGs, **R** provide the **bnlearn** package (**B**ayesian **n**etwork **learn**ing). But before using it, we have implemented the **K2** - algorithm and tested it on 3 different data sets: `Ruiz, Asia` and `Child` data sets.

# $\texttt{K2}$ algorithm

When dealing with Bayesian networks, a natural problem that arises is the determination of the most probable structure given a database $D$. Indicating by $B_S$ a network structure, we know that the probability $P(B_S|D)$ of having a particular structure given the dataset obeys

$$
\begin{equation}
P(B_S|D) = \frac{P(B_S, D)}{P(D)}
\end{equation}
$$

The challenging part is the computation of the denominator $P(D)$, since it would be necessary to evaluate all the possible network structures $\lbrace B_{S_j}\rbrace _{j \in \mathbf{N}}$. However, in practice, we are interested in ranking different structures by computing the ratio of their posterior probabilities:

$$
\begin{equation}
\frac{P(B_{S_i}|D)}{P(B_{S_j}|D)} = \frac{P(B_{S_i}, D)}{P(B_{S_j}, D)}
\end{equation}
$$

Thus, the problem of maximizing $P(B_S|D)$ can be shifted to the problem of maximizing $P(B_S,D)$. The derivation of $P(B_S,D)$ relies on the following assumptions:

1. The dataset variables are discrete.
2. Cases occur independently.
3. There are no cases that have variables with missing values.
4. Let $B_P$ denote an assignment of numerical probability values to a belief network that has structure $B_S$, the term $f(B_P|B_S)$ denotes the **likelihood** of the particular numerical probability and it is uniform. The term $P(B_S)$ can be viewed as a prior probability - namely the *preference bias probability*.

Application of Assumption 1 allows us to write 

$$
\begin{equation} P(B_S, D) = \int_{B_P} P(D|B_P, B_S) f(B_P|B_{S})P(B_S) dB_P \end{equation}
$$

with $P(D|B_P, B_S)$ being the probability mass function of our data given the network structure and its conditional probabilities assignments. Assumption 2 gives

$$
\begin{equation} P(B_S, D) = \int_{B_P}\biggl [\displaystyle\prod_{h=1}^m P(C_h|B_S, B_P) \biggl ]f(B_P|B_{S})P(B_S) dB_P \end{equation}
$$

where $m$ is the number of cases in $D$ and $C_h$ is the $h^{th}$ case in $D$. Eventually it is possible to demonstrate that Assumptions 3 and 4 lead to obtaining

$$
\begin{equation} P(B_S,D) = P(B_S)\displaystyle\prod_{i=1}^n \displaystyle\prod_{j=1}^q \frac{(r_i-1)!}{(N_{ij}+r_i-1)!} \displaystyle\prod_{k=1}^{r_i} \alpha_{ijk}!  \end{equation}
$$

The maximization of Equation (5) involves the enumeration of all possible $B_S$, which is the problem we wanted to avoid from the beginning, since it scales exponentially. However, if we assume also that **there exists an ordering** of the variables stored in $D$ and that $P(B_S)$ **is uniform**, the maximization problem relies only on finding the best configuration of the parents compatible with the given ordering:

$$
\begin{equation} \max[P(B_S,D)] = c \displaystyle\prod_{i=1}^n \max \Bigg[\displaystyle\prod_{j=1}^q \frac{(r_i-1)!}{(N_{ij}+r_i-1)!} \displaystyle\prod_{k=1}^{r_i} \alpha_{ijk}! \Bigg] \end{equation}
$$

Indeed, the whole idea of the **K2** algorithm is to heuristically perform this maximization. It is convenient to define the following function

$$
\begin{equation} f(i, \pi_i) = \displaystyle\prod_{j=1}^{q_i} \frac{(r_i-1)!}{(N_{ij}+r_i-1)!} \displaystyle\prod_{k=1}^{r_i} \alpha_{ijk}!  \end{equation}
$$

which we will refer to as the `parent_eval` function. K2 algorithm uses a greedy-search method that begins by making the assumption that a node has no parents, and then adds incrementally that parent whose addition most increases the probability of the resulting structure. When the addition of no single parent can increase the probability, we stop adding parents to the node.

<!-- The idea is to apply the above equations iteratively for every possible $B_S$. -->

* $\pi_i$ is the set of parents of node $x_i$
* $q_i=|\phi_i|$, where $\phi_i$ is the list of possible instantiations of the parents of $x_i$ in database $D$
* $r_i=|V_i|$ where $V_i$ is the list of all possible values of the attribute $x_i$
* $\alpha_{ijk}$ is number of cases in $D$ in which the attribute $x_i$ is instantiated with its $k^{th}$ value, and the parents of $x_i$ in $\pi_i$ are instantiated with the $j^{th}$ instantiation in $\phi_i$
* $N_{ij}=\displaystyle\sum_{k=1}^{r_i}\alpha_{ijk}$, so the number of instances in the database in which the parents of $x_i$ in $\pi_i$ are instantiated with the $j^{th}$ instantiation in $\phi_i$.

```r
log.parent_eval = function(i, parents = NA, df, carray) {
    # i is the place of the variable in the assigned order, i.e. the column order in the database
    R = unique(df[, i])
    r = length(R) # number of rows
    df_names = names(df)
    if (all(is.na(parents))) {
        q = 1
        N = as.data.frame(rbind(table(df[, i])))
        mat = N
    } else {
        N1 = df %>% group_by(df[parents]) %>% mutate(index=cur_group_id()) %>% group_by(index) %>% count(df[df_names[i]])
        mat = matrix(0, nrow = max(N1$index), ncol = carray[i])
        mat[cbind(data.matrix(N1[c("index", df_names[i])]))] = N1$n
    }

    logprod = rowSums(lfactorial(mat))
    ssum = rowSums(mat)
    logout = sum(lfactorial(r - 1) - lfactorial(ssum + r - 1) + logprod)
    return(logout)
}
```

This function can be computed in $\mathcal{O}(m\hspace{0.1cm}u\hspace{0.1cm}r)$ time where $u$ is the maximum number of parents that any node is permitted to have, as designated by the user and $m$ is the number of cases in the data set $D$.

For *run-time* savings results, the logarithmic version of equation (7) has been implemented, since it requires only additions and subtractions, rather than multiplications and divisions.

The whole $\texttt{K2}$ plays around the `parent_eval` function, defined in **R** as:

```r
K2_algorithm = function(n, u, D, time.info = FALSE) {
    
    # ============================================================= #
    
    # Input:
    # n = set of n ordered nodes
    # u = an upper bound on the number of parents a node may have
    # D = database containing m cases
    # time.info = if TRUE, it calculates the computation time
    
    # ============================================================= #
    
    # Output: a list with
    # A printout of the parents of the node, for each node
    # The network score evaluated with bnlearn
    
    # ============================================================= #

    if (time.info) {start = Sys.time()}

    carray = create_count_array(D)

    # Useful variables for scoring
    nodes = names(D)                            # node names
    network.dag = empty.graph(nodes=nodes)      # network DAG (Directed Acyclic Graph)

    parents = c(rep(NA, n), vector("list"))
    for (i in 1:n) {
        if (i == 1) {next}

        len_parents = 0
        P_old = log.parent_eval(i, parents[[i]], D, carray)

        OTP = TRUE # Ok To Proceed
        while (OTP & len_parents < u) {

            # let z be the node in Pred(x_i) - π_i that maximizes f(i, π_i ∪ {z});
            if (all(is.na(parents[[i]]))) {
                indexes = 1:(i-1)
            } else {
                foo = 1:(i-1)
                indexes = foo[-parents[[i]]]
            }

            if (length(indexes) == 0) {
                OTP = FALSE
                next
            }
            P = vector('numeric', length = length(indexes))
            for (t in 1:length(indexes)) {
                cand_parents = append(parents[[i]], indexes[t])
                cand_parents = cand_parents[!is.na(cand_parents)]
                P[t] = log.parent_eval(i, cand_parents, D, carray)
            }
            if (max(P) > P_old) {
                P_old = max(P)
                parents[[i]] = append(parents[[i]], indexes[which.max(P)])
                parents[[i]] = parents[[i]][!is.na(parents[[i]])]
                len_parents = length(parents[[i]])
                
                network.dag = set.arc(network.dag, from=nodes[indexes[which.max(P)]], to=nodes[i])
            } 
            else {
                OTP = FALSE
            }
        }
    }

    if (time.info) {
        end = Sys.time()
        time = difftime(end, start, units='secs')
        if (time > 60){
            cat('\nTotal execution time:', difftime(end, start, units='mins'), 'min')
        }
        else{
            cat('\nTotal execution time:', time, 's')
        }
    }

    # Network score
    network.score = score(network.dag, D)
    cat("The Network score is", network.score, "\n  ")

    return(list(parents=parents, score=network.score))
}
```

This algorithm is evaluated for $n$ times and, in the end, only the DAG with the best score is kept and saved as the best Bayesian Network predicted.

```r
K2 = function(n, u, D, seed = 12345, num.iterations = 1){
    start = Sys.time()
    set.seed(seed)

    # scoring function
    best.score = -Inf

    for(i in 1:num.iterations){

        nodes.order = if (i == 1) names(D) else sample(names(D))
        cat('order =', nodes.order, "\n")
        for(u.single in 1:(u)){

            cat('Running iteration #', i, 'u =', u.single, "\n")
            result = K2_algorithm(n, u.single, D[, nodes.order])

            if (result$score > best.score) {
                best.score = result$score
                best.order = nodes.order
                best.dag = result$parents
                best.u = u.single
            }
        }
    }
    cat(' DONE \n')

    end <- Sys.time()
    cat('\nTotal execution time:', difftime(end, start, units='mins'), 'mins\n')

    return(list(dag=best.dag, score=best.score, order=best.order, u=best.u))
}
```

Furthermore, a `get_dag` function as been implemented: it creates a `string` that, thanks to the R package **Graphviz**, allows to get a graphical representation of the Bayesian network DAG.

### Data set: `Ruiz`

The `Ruiz` data set is a very simple and immediate data set, used just for testing: the results obtained are in perfect agreement with **[3]** and the final DAG is the one expected 

$$\begin{equation*}
\bold{x_1\to x_2\to x_3}
\end{equation*}
$$

 with a computation time of $\sim 10^{-6}\hspace{0.1cm}s$ and a final score of $-20.2$.

### Data set: `Asia`

`Asia` is a randomly generated data set of the Asia Bayesian Network, an **R** data set with $10^4$ items, no missing data, and no imputation needed, ideal for **BN** computions.

In order to select the best possible architecture, we have implemented an interation function that loops through $n\gg 1$ possible combinations and select the one with the **best score**.

The execution time is of $\sim 0.9\hspace{0.1cm}s$ and the score is $-2.23\cdot 10^4$. The best architecture found is the one shown in the picture below:

<div align="center">
    <img src=images/Asia.png width=500 height=500>
</div>

The `bnlearn` library allows also to compute conditional probabilities based both on $D$ and the DAG found. For example, looking this survey, we could interested in understand the relationship between *smoking* and *bronchitis*. With `bn.fit` we can calculate the posterior probabilities of all the involved variables and easily get:

<div align="center">
    <img src=images/bronc_smoke.png width=400 height=400>
    <img src=images/xray_either.png width=400 height=400>
</div>

<!-- Moreover, we have perfomed an analysis based on **C**onditional **P**robability **Q**uery: 

<div align="center">
    <img src=images/cpq.png width=350 height=350>
</div> -->

### Data set: `Child`

The `Child` dataset contains $5\cdot 10^3$ randomly generated items with missing data (no latent variables) of the **Child Bayesian Network**. Imputation is performed, so both raw and imputed data are present.

In order to select the best possible architecture, we have implemented an interation function that loops through $n\gg 1$ possible combinations and select the one with the **best score**.

The execution time is of $\sim 6\hspace{0.1cm}s$ and the score is $-5.99\cdot10^4$. The best architecture found is the one shown in the picture below:

<div align="center">
    <img src=images/Child.png width=600 height=600>
</div>

With `bn.fit` we can calculate the posterior probabilities of all the involved variables and easily get:

<div align="center">
    <img src=images/lungparench_disease.png width=400 height=400>
    <img src=images/lungs_chestxray.png width=400 height=400>
</div>

## `bnstruct`: data set comparison

Furthermore, we have analyzed the three data sets with the **R** library $\texttt{bnstruct}$ **[4]** and compared its performances in terms of time efficiency and the final **DAG** score obtained. We see how $\texttt{bnstruct}$ is very computationally efficient, allowing to solve for the best Bayesian Network in few seconds.

We have then tested the best DAG found by the **K2** algorithm for each set and compared its **score** with the best BN structure found by $\texttt{bnstruct}$. We report the results in the table below:

<center>

| Data set | $\tau_{iter}^\texttt{K2}$ | $score_{\texttt{K2}}$ | $score_{\texttt{bnstruct}}$ |
| :-------: | :------------------: | :---------------------: | :---------------------------: |
| `Ruiz` |    $\sim 10^{-6}\hspace{0.1cm}s$    |         $-20.2$         |               $-20.2$               |
| `Asia` |    $\sim 0.9\hspace{0.1cm}s$    |   $-2.23\cdot 10^4$   |               $-2.23\cdot 10^4$            |
| `Child` |    $\sim 6\hspace{0.1cm}s$    |    $-5.99\cdot10^4$   |               $-6.10\cdot10^4$               |

</center>

The score found by $\texttt{bnstruct}$ is different from the one found by $\texttt{K2}$.
The **R** library is able to find another DAG structure, both on small and larger data sets: this is due to the `mmhc` algorithm that performs a statistical sieving of the search space followed by a greedy evaluation. The Min-Max Hill Climb algorithm provides a considerably fast exectution time at the expense of a lower quality result **[5]**.  

The increasing execution time can be explained by considering that the $\texttt{K2}$ algortihm as an overall time complexity of $\mathcal{O}(m + r - 1) + \mathcal{O}(m u n r) \mathcal{O}(u) n = \mathcal{O}(m r u^2 n^2)$. In the worst case, when $u = n$, the time complexity is $\mathcal{O}(m r n^4)$.

Equation (1) represents more accurate and fundamental results **[2]** than the K2 algorithm does, but nonetheless K2 allows to get a preliminar grasp on the structure of the network under study.

## What's next?!

A further development than can be performed could be to implment a sort of **backward K2** algorithm in which, basically, instead of adding parents, we start with a *fully-connected* graph and interatively remove links according to the *score* function. However, K2 and the backward K2 algorithms are not mutually exclusive, but they could be run together to formulate different many structures and

## **References**

**[1]** M. Scutari and J. B. Denis, *Bayesian Networks*, CRC Press, 2022, Taylor and Francis Group,

**[2]** G. F. Cooper and E. Herskovits, *A Bayesian Method for the Induction of Probabilistic Networks from Data*, Machine Learning 9, (1992) 309,

**[3]** C. Ruiz, *Illustration of the K2 Algorithm for learning Bayes Net Structures*, http://web.cs.wpi.edu/~cs539/s11/Projects/k2_algorithm.pdf",

**[4]** A. Franzin et al., *$\texttt{bnstruct}$: an R package for Bayesian Network structure learning in the presence of missing data*, Bioinformatics 33(8) (2017) 1250,

**[5]** F. Sambo and A. Franzin, *$\texttt{bnstruct}$: an R package for Bayesian Network Structure Learning with missing data*, December 12, 2016,
