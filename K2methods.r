library(bnlearn)

######################################################################################
# parent_eval
######################################################################################
# calcualtes function (12) in the article
#
#
#
#
#
#
#
######################################################################################
parent_eval = function(i, parents = NA, df) { # i is the place of the variable in the assigned order, not the order of the variable in the database
    R = unique(df[, i])
    r = length(R)

    if (all(is.na(parents))) {
        q = 1
        N = as.data.frame(rbind(table(df[, i])))
    } else {
        W = as.data.frame(cbind(unique(df[, parents])))
        q = nrow(W)
        N = data.frame(matrix(nrow = q, ncol = r))

        for (j in 1:q) {
            df_parents = df[apply(as.data.frame(cbind(df[,parents])), 1, function(row) all(row == as.numeric(W[j,]))),]
            for (k in 1:r) {
                N[j, k] = sum(df_parents[, i] == R[k])
            }
        }
    }

    # cat('matrix N referred to variable i:', i, '\n')
    # print(N)

    # here the lfactorial could give problems with small datasets
    foo = data.frame(prod = apply(as.data.frame(map(N, lfactorial)), 1, prod), 
                    sum = apply(N, 1, sum))

    # cat('\ndataframe foo referred to variable i:', i, '\n')
    # print(foo)

    # here the lfactorial at the denominator could give problems with small datasets
    out = prod(unlist(map2(.x = foo$prod, .y = foo$sum, .f = ~ (factorial(r - 1) / lfactorial(.y + r - 1)) * .x)))
    return(out)
}

######################################################################################
# K2_algorithm
######################################################################################
# implements K2 algorithm
#
#
#
#
#
#
#
######################################################################################
K2_algorithm = function(n, u, D) {
    parents = c(rep(NA, n), vector("list"))
    for (i in 1:n) {
        if (i == 1) {next}

        len_parents = 0
        P_old = parent_eval(i, parents[[i]], D)
        # cat('\nPold:', P_old, 'i:', i, '\n')

        OTP = TRUE
        while(OTP & len_parents < u) {
            # cat('parents:', parents[[i]], 'i:', i, '\n')
            # cat('condition:', is.na(parents[[i]]), 'variable:', i, '\n')
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
                P[t] = parent_eval(i, cand_parents, D)
                # cat('P[t]:', P[t], 't:', t, 'i:', i, '\n')
            }
            # cat('P vector:', P, 'i:', i, '\n')
            # cat('maxP:', max(P), 'Pold:', P_old, 'i:', i, 'parents:', parents[[i]], '\n')
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
    return(parents)
}

######################################################################################
# get_dag
######################################################################################
# returns an object handable by the class bnlearn
#
#
#
#
#
#
#
######################################################################################
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
    # cat('the string is:', string, '\n')
    return(model2network(string))
}