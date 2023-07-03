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

Consider the problem of determining a belief-network structure $B_S$ that maximizes
$P(B_S|D)$. In general, more than a single structure can be reasonable, but we will assume here that the structure that maximizes the probability is only one. For a given data set $D$, $P(B_S, D) \propto P(B_S|D)$, so that finding the $B_S$ that maximizes $P(B_S|D)$ is equivalent to find the $B_S$ that maximizes $P(B_S, D)$.

**Assumptions:**

1. The data set variables are discrete.
2. Cases occur independently.
3. There are no cases that have variables with missing values.
4. Let $B_P$ denote an assignment of numerical probability values to a belief network that has structure $B_{S_1}$, the term $f(B_P|B_{S_1})$ denotes the **likelihood** of the of the particular numerical probability and it is uniform. The term $P(B_S)$ can be viewed as a prior probability - namely the *preference bias probability*.

The whole idea of the **K2** algorithm is to heuristically look for the most probable belief–network structure given a database of cases. This is done, by computing 

$$
\begin{equation} P(B_S, D) = P(B_S) \displaystyle\prod_{i=1}^n f(i, \pi_i) \end{equation}
$$
where $f(B_P|B_S)$ can be expressed (**Theorem 1**) as 

$$
\begin{equation} f(i, \pi_i) = \displaystyle\prod_{j=1}^q \frac{(r_i-1)!}{(N_{ij}+r_i-1)!} \displaystyle\prod_{k=1}^{r_i} \alpha_{ijk}!  \end{equation}
$$

 and we refer to it as the `parent_eval` function (computed as logarithm in order to work with sums rather than products). This is a heuristic method that uses a greedy-search method that begins by making the assumption that a node has no parents, and then adds incrementally that parent whose addition most increases the probability of the resulting structure. When the addition of no single parent can increase the probability, we stop adding parents to the node.

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
        mat[cbind( data.matrix(N1[c("index", df_names[i])]))] = N1$n

    }

    logprod = rowSums(lfactorial(mat))
    ssum = rowSums(mat)
    logout2 = sum(lfactorial(r - 1) - lfactorial(ssum + r - 1) + logprod)
    return(logout2)
}
```

This function can be computed in $\mathcal{O}(m\:u\:r)$ time where $u$ is the maximum number of parents that any node is permitted to have, as designated by the user and $m$ is the number of cases in the data set $D$.

For *run-time* savings results, the logarithmic version of equation (2) has been implemented, since it requires only additions and subtractions, rather than multiplications and divisions.

The whole $\texttt{K2}$ plays around the `parent_eval` function: we have defined it in **R** as:

```r
K2_algorithm = function(n, u, D, time.info = FALSE) {
    
    # ============================================================= #
    
    # Input:
    # n = set of n ordered nodes
    # u = an upper bound on the number of parents a node may have
    # D = database containing m cases
    # time.info = if TRUE, it calculates the computation time
    
    # ============================================================= #
    
    # Output: 
    # For each node, a printout of the parents of the node
    
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
        # cat('\nPold:', P_old, 'i:', i, '\n')

        OTP = TRUE # Ok To Proceed
        while (OTP & len_parents < u) {

            # let z be the node in Pred(x_i) - π_i that maximizes f(i, π_i ∪ {z});
            if (all(is.na(parents[[i]]))) {
                indexes = 1:(i-1)
            } else {
                foo = 1:(i-1)
                indexes = foo[-parents[[i]]]
            }
            # cat('indexes:', indexes, 'i:', i, '\n')

            if (length(indexes) == 0) {
                OTP = FALSE
                next
            }
            P = vector('numeric', length = length(indexes))
            for (t in 1:length(indexes)) {
                cand_parents = append(parents[[i]], indexes[t])
                cand_parents = cand_parents[!is.na(cand_parents)]
                # cat('cand_parents:', cand_parents, 't:', t, 'i:', i, '\n')
                P[t] = log.parent_eval(i, cand_parents, D, carray)
                # cat('P[t]:', P[t], 't:', t, 'i:', i, '\n')
            }
            # cat('P vector:', P, 'i:', i, '\n')
            # cat('maxP:', max(P), 'Pold:', P_old, 'i:', i, 'parents:', parents[[i]], '\n')
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

The `Ruiz` data set is a very sample and immediate data set, used just for testing: the results obtained are in perfect agreement with **[3]** and the final DAG is the one expected 

$$
\bold{x_1\to x_2\to x_3}
$$

 with a computation time of $\sim 0.006\:s$ and a final score of $-20$.

### Data set: `Asia`

`Asia` is a randomly generated data set of the Asia Bayesian Network, an **R** data set with $10^4$ items, no missing data, and no imputation needed, ideal for **BN** computions.

In order to select the best possible architecture, we have implemented an interation function that loops through *n>1* possible combinations and select the one with the **best score**.

The execution time is of $\sim 30\:mins$ and the score is $-2.2\cdot 10^4$. The best architecture found is the one shown in the picture below:

<div align="center">
    <img src=images/Asia.png width=400 height=400>
</div>

### Data set: `Child`

The `Child` dataset contains $5\cdot 10^3$ randomly generated items with missing data (no latent variables) of the **Child Bayesian Network**. Imputation is performed, so both raw and imputed data are present.

In order to select the best possible architecture, we have implemented an interation function that loops through *n=1* possible combinations and select the one with the **best score**.

The execution time is of $\sim 70\:mins$ and the score is $-6.0\cdot10^4$. The best architecture found is the one shown in the picture below:

<div align="center">
    <img src=images/Child.png width=450 height=450>
</div>

## `bnstruct` data set comparison

Furthermore, we have analyzed the three data sets with the **R** library $\texttt{bnstruct}$ **[4]** and compared its performances in terms of time efficiency and the final **DAG** score obtained. We see how $\texttt{bnstruct}$ is very computationally efficient, allowing to solve for the best Bayesian Network in few seconds.

We have then tested the best DAG found by the **K2** algorithm for each set and compared its **score** with the best BN structure found by $\texttt{bnstruct}$. We report the results in the table below:

| Data set | $\tau_{iter}^{K2}$ | $score_{\texttt{K2}}$ | $score_{\texttt{bnstruct}}$ |
| :-------: | :------------------: | :---------------------: | :---------------------------: |
| `Ruiz` |    $\sim 1\:s$    |         $-0.2\cdot 10^2$         |               ?               |
| `Asia` |    $\sim 5\:m$    |   $-2.2\cdot 10^4$   |               ?               |
| `Child` |    $\sim 1\:h$    |    $-6.0\cdot10^4$    |               ?               |

The score found $\texttt{bnstruct}$ is different from the one found by $\texttt{K2}$: this means that the **R** library is able to find a better DAG structure, both on small and larger dataset.

The very huge exectution time can be explained by considering that the $\texttt{K2}$ algortihm as an overall time complexity of $\mathcal{O}(m + r - 1) + \mathcal{O}(m u n r) \mathcal{O}(u) n = \mathcal{O}(m r u^2 n^2)$. In the worst case, when $u = n$, the time complexity is $\mathcal{O}(m r n^4)$.

Equation (1) represents more accurate and fundamental results **[2]** than the K2 algorithm does, but nonetheless K2 allows to get a preliminar grasp on the structure of the network under study.

## What's next?!

A further development than can be performed could be to implment a sort of **backward K2** algorithm in which, basically, instead of adding parents, we start with a *fully-connected* graph and interatively remove links according to the *score* function. However, K2 and the backward K2 algorithms are not mutually exclusive, but they could be run together to formulate different many structures and

## **References**

**[1]** M. Scutari and J. B. Denis, *Bayesian Networks*, CRC Press, 2022, Taylor and Francis Group,

**[2]** G. F. Cooper and E. Herskovits, *A Bayesian Method for the Induction of Probabilistic Networks from Data*, Machine Learning 9, (1992) 309,

**[3]** C. Ruiz, *Illustration of the K2 Algorithm for learning Bayes Net Structures*, http://web.cs.wpi.edu/~cs539/s11/Projects/k2_algorithm.pdf",

**[4]** A. Franzin et al., *$\texttt{bnstruct}$: an R package for Bayesian Network structure learning in the presence of missing data*, Bioinformatics 33(8) (2017) 1250,

**[5]** F. Sambo and A. Franzin, *$\texttt{bnstruct}$: an R package for Bayesian Network Structure Learning with missing data*, December 12, 2016,
