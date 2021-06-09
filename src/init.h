#ifndef PHYLO_INIT_H
#define PHYLO_INIT_H

#include <R.h>
#include <Rinternals.h>

/* treeio.c */
SEXP phylo_phy_read_newickstr(SEXP);
SEXP phylo_phy_write_newickstr(SEXP);
SEXP phylo_tiplabels(SEXP);
SEXP phylo_node_notes(SEXP);
SEXP phylo_phy_node_brlens(SEXP);
SEXP phylo_phy_node_ages(SEXP);
SEXP phylo_phy_node_ancestors(SEXP, SEXP);
SEXP phylo_phy_node_children(SEXP, SEXP);
SEXP phylo_phy_node_descendants(SEXP, SEXP, SEXP, SEXP);
SEXP phylo_phy_extract_clade(SEXP, SEXP);
SEXP phylo_phy_extract_subtree(SEXP, SEXP, SEXP);
SEXP phylo_phy_ladderize(SEXP);
SEXP phylo_phy_node_rotate(SEXP, SEXP);
/* treeplot.c */
SEXP phylo_plot_cartesian(SEXP, SEXP, SEXP);
SEXP phylo_plot_polar(SEXP, SEXP);

#endif
