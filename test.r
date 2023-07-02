# source("K2methods.r")
suppressWarnings(suppressMessages(library(purrr)))
suppressWarnings(suppressMessages(library(Rgraphviz)))
suppressWarnings(suppressMessages(library(bnlearn)))
suppressWarnings(suppressMessages(library(bnstruct)))
suppressWarnings(suppressMessages(library(dplyr)))
# set.seed(9112001)

# data <- NULL
# data$case <- as.factor(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)+1)
# data$x1 <- as.factor(c(1, 1, 0, 1, 0, 0, 1, 0, 1, 0)+1)# classification target
# data$x2 <- as.factor(c(0, 1, 0, 1, 0, 1, 1, 0, 1, 0)+1)
# data$x3 <- as.factor(c(0, 1, 1, 1, 0, 1, 1, 0, 1, 0)+1)
# data <- as.data.frame(data, stringsAsFactors = TRUE)

# data = asia()
# data = impute(data)
# data = imputed.data(data)
#
# mode(data) = 'character'
# DF = data.frame(data, stringsAsFactors = TRUE)
# cat('The dataset dimensions are:', dim(DF))
#
# colnames(DF) <- asia()@variables
# head(DF)
# today = 01072023 # default seed = 1234
# struct = K2(n = ncol(DF), u = 3, D = DF, seed=today, num.iterations=1)

create_count_array = function(df) {
    df_names = names(df)
    array = array(0, dim = ncol(df))
    for (i in 1:length(df_names)) {
        array[i] = nrow(unique(df[df_names[i]]))
    }
    return(array)
}


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
            cat('\nTotal execution time:', time, 's', "\n")
        }
    }

    # Network score
    network.score = score(network.dag, D)
    # cat("The Network score is", network.score, "\n  ")

    return(list(parents=parents, score=network.score))
}


K2 = function(n, u, D, seed = 12345, num.iterations = 1){
    start = Sys.time()
    set.seed(seed)

    # scoring function
    best.score = -Inf

    for(i in 1:num.iterations){

        nodes.order = if (i == 1) names(D) else sample(names(D))
        cat('order =', nodes.order)

        cat('Running iteration #', i, 'u =', u, "\n")
        result = K2_algorithm(n, u, D[, nodes.order])

        if (result$score > best.score) {
            best.score = result$score
            best.order = nodes.order
            best.dag = result$parents
            best.u = u
        }
    }
    cat(' DONE \n')

    end <- Sys.time()
    cat('\nTotal execution time:', difftime(end, start, units='mins'), 'mins\n')

    return(list(dag=best.dag, score=best.score, order=best.order, u=best.u))
}

get_gigi_dag = function(vars, parents) {
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

# main = function(){
#     data = asia()
#     data = impute(data)
#     data = imputed.data(data)
#
#     # mode(data) = 'character'
#     DF = data.frame(data, stringsAsFactors = TRUE)
#     cat('The dataset dimensions are:', dim(DF))
#
#     colnames(DF) <- asia()@variables
#     head(DF)
#
#     K2_algorithm(n = ncol(DF), u = 3, D = DF)
# }
