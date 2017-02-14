#ifndef _OPTIONS_H
#define _OPTIONS_H

typedef struct options_t{
  int verbose;
  int ignore;
  int norm_order;
  int reorder;
  double mesh_tol;
  double err_tol;
  double perr_tol;
  char *file1, *file2;
  char **extra_varnames;
  int n_extra_varnames;
  char **ignorableVarnames;
  int n_ignorableVarnames;
  int gridVarSelfDiscovery;
  char benchmark;
}options_t;

int parse_cmdline(options_t *opts, int argc, char **argv);

#endif /* _OPTIONS_H */
