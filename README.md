<!-- # R - K2 algorithm -->

<h1 align="center">Advanced Statistics - PoD<br>AY 2022/2023 <br> University of Padua</h1>

<p align="center">
  <img src="https://user-images.githubusercontent.com/62724611/166108149-7629a341-bbca-4a3e-8195-67f469a0cc08.png" height="150"/>
  &emsp;
  <img src="https://user-images.githubusercontent.com/62724611/166108076-98afe0b7-802c-4970-a2d5-bbb997da759c.png" height="150"/>
</p>

<h3 align="center"><b>Group members:</b> E. Sarte, R. Tancredi, G. Zago<br></h3>

# Bayesian-Networks
A **Bayesian Network**, known also as *belief network*, is a probabilistic graphical model that represents the relationships between given variables as a *directed acyclic graph* - **DAG**s.

Each **node** in the graph corresponds to one of the variables $A, B, C, D \dots$ in the survey and the *direct* and *indirect* relations are expressed via **arcs** between pairs of variables: $A\to E$ means that $E$, the *child* node, depends on $A$, the *parent*.

To create and manipulating BN DAGs, **R** provide the **bnlearn** package (**B**ayesian **n**etwork **learn**ing). But before using it, we have implemented the **K2** - algorithm and tested it on 3 different data sets: `Ruiz, Asia` and `Child` data sets.

# K2 algorithm

Consider the problem of determining a belief-network structure $B_S$ that maximizes
$P(B_S|D)$. In general, more than a single structure can be reasonable, but we will assume here that the structure that maximizes the probability is only one. For a given data set $D$, $P(B_S, D) \propto P(B_S|D)$, so that finding the $B_S$ that maximizes $P(B_S|D)$ is equivalent to find the $B_S$ that maximizes $P(B_S, D)$. 

**Assumptions:**
1. The data set variables are discrete.
2. Cases occur independently.
3. There are no cases that have variables with missing values.
4. Let $B_P$ denote an assignment of numerical probability values to a belief network that has structure $B_{S_1}$, the term $f(B_P|B_{S_1})$ denotes the **likelihood** of the of the particular numerical probability and it is uniform. The term $P(B_S)$ can be viewed as a prior probability - namely the *preference bias probability*. 

The whole idea of the **K2** algorithm is to heuristically look for the most probable belief–network structure given a database of cases. This is done, by computing $$ P(B_S, D) = P(B_S) \displaystyle\prod_{i=1}^n f(i, \pi_i) $$ where $f(B_P|B_S)$ can be expressed (**Theorem 1**) as $$ f(i, \pi_i) = \displaystyle\prod_{j=1}^q \frac{(r_i-1)!}{(N_{ij}+r_i-1)!} \displaystyle\prod_{k=1}^{r_i} \alpha_{ijk}! $$ and we refer to it as the `parent-eval` function (computed as logarithm in order to work with sums rather than products). This is a heuristic method that uses a greedy-search method that begins by making the assumption that a node has no parents, and then adds incrementally that parent whose addition most increases the probability of the resulting structure. When the addition of no single parent can increase the probability, we stop adding parents to the node. 
<!-- The idea is to apply the above equations iteratively for every possible $B_S$. -->
* $\pi_i$ is the set of parents of node $x_i$
* $q_i=|\phi_i|$, where $\phi_i$ is the list of possible instantiations of the parents of $x_i$ in database $D$
* $r_i=|V_i|$ where $V_i$ is the list of all possible values of the attribute $x_i$
* $\alpha_{ijk}$ is number of cases in $D$ in which the attribute $x_i$ is instantiated with its $k^{th}$ value, and the parents of $x_i$ in $\pi_i$ are instantiated with the $j^{th}$ instantiation in $\phi_i$
* $N_{ij}=\displaystyle\sum_{k=1}^{r_i}\alpha_{ijk}$, so the number of instances in the database in which the parents of $x_i$ in $\pi_i$ are instantiated with the $j^{th}$ instantiation in $\phi_i$.

```r
log.parent_eval = function(i, parents = NA, df) { 
    # i is the place of the variable in the assigned order, i.e. the column order in the database
    R = unique(df[, i])
    r = length(R) # Number of rows

    if (all(is.na(parents))) {
        q = 1
        N = as.data.frame(rbind(table(df[, i])))
    } 
    else {
        W = unique(df[, parents, drop = FALSE])
        q = nrow(W)
        N = data.frame(matrix(nrow = q, ncol = r))
        for (j in 1:q) {
            df_parents = df[apply(df[, parents, drop = FALSE], 1, function(row) all(row == W[j, ])), ]
            N[j, ] = as.data.frame(rbind(table(df_parents[, i])))   # this automatically sorts the counts per variables per instantiation 
                                                                    # by increasing variable value (i.e. 1 count for variable 0, 1 count
                                                                    # for variable 1 etc etc)
        }
    }
    # here the lfactorial could give problems with small datasets:
    # to avoid it we compute the log(f) function:                            
    log.foo = data.frame(logprod = pmap_dbl(map_df(N, lfactorial), sum), 
                    sum = pmap_dbl(N, sum))

    log.out = sum(map2_dbl(.x = log.foo$logprod, .y = log.foo$sum, .f = ~ (lfactorial(r - 1) - lfactorial(.y + r - 1)) + sum(.x)) )
    return(log.out)
}
```

This function can be computed in $\mathcal{O}(m\:u\:r)$ time where $u$ is the maximum number of parents that any node is permitted to have, as designated by the user and $m$ is the number of cases in the data set $D$.  

This function is the essence of the algorithm, which is instantiated in **R** as: 

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
    parents = c(rep(NA, n), vector("list"))
    for (i in 1:n) {
        if (i == 1) {next}

        len_parents = 0
        P_old = log.parent_eval(i, parents[[i]], D)

        OTP = TRUE # Ok To Proceed
        while (OTP && len_parents < u) {

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
                P[t] = log.parent_eval(i, cand_parents, D)
            }
            if (max(P) > P_old) {
                P_old = max(P)
                parents[[i]] = append(parents[[i]], indexes[which.max(P)])
                parents[[i]] = parents[[i]][!is.na(parents[[i]])]
                len_parents = length(parents[[i]])
            } else {
                OTP = FALSE
            }
        }
    }

    if (time.info) {
        end = Sys.time()
        cat('\nTotal execution time:', difftime(end, start, units='secs'), 's')    
    }

    return(parents)
}
```

Furthermore, a `get_dag` function as been implemented: it creates a `string` that,thanks to the R package **Graphviz**, allows to get a graphical representation of the Bayesian network DAG.

### Data set: `Ruiz`

The `Ruiz` data set is a very sample and immediate data set, used just for testing: the results obtained are in perfect agreement with **[3]** and the final DAG is the one expected $$\bold{x_1\to x_2\to x_3}$$ with a computation time of $\sim 0.006\:s$ and a final score of $-20$. 

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

## **References**
**[1]** M. Scutari and J. B. Denis, *Bayesian Networks*, CRC Press, 2022, Taylor and Francis Group,

**[2]** G. F. Cooper and E. Herskovits, *A Bayesian Method for the Induction of Probabilistic Networks from Data*, Machine Learning 9, (1992) 309,

**[3]** C. Ruiz, *Illustration of the K2 Algorithm for learning Bayes Net Structures*, http://web.cs.wpi.edu/~cs539/s11/Projects/k2_algorithm.pdf",

**[4]** A. Franzin et al., *$\texttt{bnstruct}$: an R package for Bayesian Network structure learning in the presence of missing data*, Bioinformatics 33(8) (2017) 1250,

**[5]** F. Sambo and A. Franzin, *$\texttt{bnstruct}$: an R package for Bayesian Network Structure Learning with missing data*, December 12, 2016,
