for(i in 1:27) {
  cmd <- paste0('sbatch -p exec -t 1200 --cpus-per-task 1 --mem 16000 --wrap="JOB_IND=',i,' Rscript ./inst/raw/create_variance_from_additive_lookup_table_step1.R"')
  print(cmd)
  system(cmd)
}
