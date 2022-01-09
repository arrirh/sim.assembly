#' Performes community modeling and scaling of phylogenetic and functional structure
#'
#' @param script A Dataframe containing parameters of model and scaling for each script.
#' @param n_iter A number of iterations for each script.
#' @param sim_verb_iter A number of iterations in community simulation after which the message is returned.
#' @param be_verbose Boolean, returns messages during the execution.
#' @param make_backup Boolean, indicates if results should be restored in backup file.
#' @param use_backup Boolean, use a previous backup after interuption.
#' @param backup_path A character, path of backup file.
#'
#' @return List of generated communities and scaling results.
#' @export
#'
#' @examples
#' scaling_pipeline(script, n_iter, sim_verb_iter = 100, be_verbose = T,
#' make_backup = T, use_backup = F, backup_path = "backup.rda")


scaling_pipeline <- function(script, n_iter, sim_verb_iter = 100, be_verbose = T,
                             make_backup = T, use_backup = F, backup_path = "backup.rda"){

  if (use_backup){
    print('Restoring from backup')

    e <- new.env()
    load(backup_path, envir = e)
    e

    ii_start <- ii
    jj_start <- jj

  } else {

    coms <- replicate(n = nrow(script), expr = list())

    empty <- data.frame(r = numeric(), ses = numeric())

    res <- list(tr = list(), phy = list())
    res$tr <- vector("list", nrow(script))
    res$phy <- vector("list", nrow(script))

    for (ii in 1:nrow(script)) {
      res$tr[[ii]] <- list(G = empty, P = empty, C = empty)
      res$phy[[ii]] <- list(G = empty, P = empty, C = empty)
    }

    jj_start <- 1
    ii_start <- 1
  }

  jj_array <- jj_start:n_iter

  for (jj in jj_array) { #позволяет начать с того ii, с которого остановились в бэкапе. если итерация jj, на котором все закончилось, прошла, то ii снова начинается с 1
    if (jj > jj_start){
      ii_start <- 1
    }
    for (ii in ii_start:nrow(script)) {

      if (be_verbose){
        print(paste("iteration", jj, "/ script", ii))
      }

      coms[[ii]][[jj]] <- create_com(arena_size = script$arena_size[ii],
                                     arena_roughness = script$arena_roughness[ii],
                                     phylo_b = script$phylo_b[ii],
                                     phylo_d = script$phylo_d[ii],
                                     phylo_n_taxa = script$phylo_n_taxa[ii],
                                     trait_model = script$trait_model[ii],
                                     trait_sigma = script$trait_sigma[ii],
                                     trait_alpha = script$trait_alpha[ii],
                                     com_range = script$com_range[ii],
                                     betaEnv = script$betaEnv[ii],
                                     betaComp = script$betaComp[ii],
                                     com_nu = script$com_nu[ii],
                                     com_sigma = script$com_sigma[ii],
                                     sim_n_iter = script$sim_n_iter[ii],
                                     be_verbose = be_verbose,
                                     sim_verb_iter = sim_verb_iter)

      temp <- com_scaling(com = coms[[ii]][[jj]],
                          grain_size = script$grain_size[ii],
                          use_tr = T,
                          abundance_weighted = script$abundance_weighted[ii],
                          scaling_iter = script$scaling_iter[ii],
                          be_verbose = be_verbose)

      res$tr[[ii]]$G <- rbind(res$tr[[ii]]$G, temp$resG)
      res$tr[[ii]]$P <- rbind(res$tr[[ii]]$P, temp$resP)
      res$tr[[ii]]$C <- rbind(res$tr[[ii]]$C, temp$resC)


      temp <- com_scaling(com = coms[[ii]][[jj]],
                          grain_size = script$grain_size[ii],
                          use_tr = F,
                          abundance_weighted = script$abundance_weighted[ii],
                          scaling_iter = script$scaling_iter[ii])

      res$phy[[ii]]$G <- rbind(res$phy[[ii]]$G, temp$resG)
      res$phy[[ii]]$P <- rbind(res$phy[[ii]]$P, temp$resP)
      res$phy[[ii]]$C <- rbind(res$phy[[ii]]$C, temp$resC)


      if (make_backup){
        save(coms, res, ii, jj, file = backup_path)
      }

    }
  }

  final_res <- list(coms = coms, res = res)

  return(final_res)

}
