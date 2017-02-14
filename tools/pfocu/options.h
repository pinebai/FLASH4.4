#ifndef _OPTIONS_H
#define _OPTIONS_H

typedef struct options_t{
  int verbose;
  int ignore;
  double mesh_tol;
  double err_tol;
  double perr_tol;
  char *file1, *file2;
  char **extra_varnames;
  int n_extra_varnames;
  char **ignorableVarnames;
  int n_ignorableVarnames;
}options_t;

int parse_cmdline(options_t *opts, int argc, char **argv);

#endif /* _OPTIONS_H */
