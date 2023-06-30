library(bnlearn)
library(data.table)
library(dplyr)

#######################
### parent_eval     ###
### formula (12)    ###
### in the article  ###
#######################

parent_eval = function(i, parents = NA, df) { # i is the place of the variable in the assigned order, i.e. the column order in the database
    R = unique(df[, i])
    r = length(R)

    if (all(is.na(parents))) {
        q = 1
        N = as.data.frame(rbind(table(df[, i])))
    } else {
        W = unique(df[, parents, drop = FALSE])
        q = nrow(W)
        N = data.frame(matrix(nrow = q, ncol = r))

        for (j in 1:q) {
            # mask = pmap_lgl(DF[, c(1,2), drop = FALSE], function(...) all(c(...) == comp_vec))
            # DF[mask,,drop = FALSE]
            df_parents = df[apply(df[, parents, drop = FALSE], 1, function(row) all(row == as.numeric(W[j, ]))), ]
            N[j, ] = as.data.frame(rbind(table(df_parents[, i]))) # this automatically sorts the counts per variables per intstantiation 
                                                                    # by increasing variable value (i.e. 1 count for variable 0, 1 count
                                                                    # for variable 1 etc etc)
        }
    }

    # cat('matrix N referred to variable i:', i, '\n')
    # print(N)

    # here the lfactorial could give problems with small datasets
    foo = data.frame(prod = pmap_dbl(map_df(N, lfactorial), prod), 
                    sum = pmap_dbl(N, sum))

    # cat('\ndataframe foo referred to variable i:', i, '\n')
    # print(foo)

    # here the lfactorial at the denominator could give problems with small datasets
    out = prod(map2_dbl(.x = foo$prod, .y = foo$sum, .f = ~ (factorial(r - 1) / lfactorial(.y + r - 1)) * .x))
    return(out)
}

#######################
### log.parent_eval ###
### formula (12)    ###
### in the article  ###
#######################
                                  
log.parent_eval = function(i, parents = NA, df) { 
    # i is the place of the variable in the assigned order, i.e. the column order in the database
    R = unique(df[, i])
    r = length(R) # number of rows
    # cat("Parents are:\n")
    # print(parents)
    if (all(is.na(parents))) {
        q = 1
        N = as.data.frame(rbind(table(df[, i])))
    } else {
        W = unique(df[, parents, drop = FALSE])
        q = nrow(W)
        N = data.frame(matrix(nrow = q, ncol = r))
        # cat("W matrix is:\n")
        # print(W)
        for (j in 1:q) {
            # mask = pmap_lgl(DF[, c(1,2), drop = FALSE], function(...) all(c(...) == comp_vec))
            # DF[mask,,drop = FALSE]
            # cat("Exploring W elemnts:\n")
            # print(W[j, , drop = TRUE])
            # cat("df[, parents]:\n")
            # print(df[, parents, drop = FALSE])
            # cat("long function:\n")
            # print(apply(df[, parents, drop = FALSE], 1, function(row) all(row == as.numeric(j-1))))
            df_parents = df[apply(df[, parents, drop = FALSE], 1, function(row) all(row == W[j, ])), ]
            # cat("df_parents is:\n")
            # print(df_parents)
            N[j, ] = as.data.frame(rbind(table(df_parents[, i])))   # this automatically sorts the counts per variables per instantiation 
                                                                    # by increasing variable value (i.e. 1 count for variable 0, 1 count
                                                                    # for variable 1 etc etc)
        }
    }

    # cat('matrix N referred to variable i:', i, '\n')
    # print(N)

    # here the lfactorial could give problems with small datasets:
    # to avoid it we compute the log(f) function:                            
    log.foo = data.frame(logprod = pmap_dbl(map_df(N, lfactorial), sum), 
                    sum = pmap_dbl(N, sum))

    # cat('\ndataframe foo referred to variable i:', i, '\n')
    # print(log.foo)

    # here the lfactorial at the denominator could give problems with small datasets
#     out = prod(map2_dbl(.x = foo$prod, .y = foo$sum, .f = ~ (factorial(r - 1) / lfactorial(.y + r - 1)) * .x))
    log.out = sum(map2_dbl(.x = log.foo$logprod, .y = log.foo$sum, .f = ~ (lfactorial(r - 1) - lfactorial(.y + r - 1)) + sum(.x)) )
    # cat(log.out)
    return(log.out)
}

####################
### K2_algorithm ###
####################

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

    # Useful variables for scoring
    nodes = names(D)                     # node names
    network.dag = empty.graph(nodes=nodes)   # network DAG (Directed Acyclic Graph)

    parents = c(rep(NA, n), vector("list"))
    for (i in 1:n) {
        if (i == 1) {next}

        len_parents = 0
        P_old = log.parent_eval(i, parents[[i]], D)
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
                P[t] = log.parent_eval(i, cand_parents, D)
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

K2 <- function(n, u, D, seed = 12345, num.iterations = 1){
    start = Sys.time()
    set.seed(seed)

    # scoring function
    best.score = -Inf

    for(i in 1:num.iterations){
        
        cat('Running iteration #', i, '...')

        nodes.order = if (i == 1) names(D) else sample(names(D))
        result = K2_algorithm(n, u, D)

        if (result$score > best.score) {
            best.score = result$score
            best.order = nodes.order
            best.dag = result$parents
        }
    }
    cat(' DONE \n')

    end <- Sys.time()
    cat('\nTotal execution time:', difftime(end, start, units='mins'), 'mins\n')

    return(list(dag=best.dag, score=best.score, order=best.order))
}

###########################
### get_dag             ###
### returns an object   ### 
### handable by the     ###
### bnlearn class       ###
###########################

get_dag = function(vars, parents) {
    string = ''
    for (i in 1:length(vars)) {
        string = paste(string, '[', vars[i], sep = '')
        if (all(is.na(parents[[i]]))) {
            string = paste(string, ']', sep = '')
        } else {
            if (length(parents[[i]]) == 1) {
                string = paste(string, '|', vars[parents[[i]]], sep = '')
            } else {
                for (j in 1:length(parents[[i]])) {
                    if (j == 1) {
                        string = paste(string, '|', vars[parents[[i]]][j], ':', sep = '')
                    } else if (j == length(parents[[i]])) {
                        string = paste(string, vars[parents[[i]]][j], sep = '')
                    } else {
                        string = paste(string, vars[parents[[i]]][j], ':', sep = '')
                    }
                }
            }

           string = paste(string, ']', sep = '')
        }
    }
#     cat('the string is:', string, '\n')
    return(model2network(string))
}