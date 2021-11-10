#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>

#include "phy.h"
#include "init.h"

#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}


static const R_CallMethodDef CallEntries[] = {
    CALLDEF(phylo_phy_read_newickstr, 1),
    CALLDEF(phylo_phy_write_newickstr, 1),
    CALLDEF(phylo_tiplabels, 1),
    CALLDEF(phylo_node_notes, 1),
    CALLDEF(phylo_phy_node_brlens, 1),
    CALLDEF(phylo_phy_node_ages, 1),
    CALLDEF(phylo_phy_node_ancestors, 2),
    CALLDEF(phylo_phy_node_children, 2),
    CALLDEF(phylo_phy_node_descendants, 4),
    CALLDEF(phylo_phy_extract_clade, 2),
    CALLDEF(phylo_phy_extract_subtree, 3),
    CALLDEF(phylo_phy_ladderize, 1),
    CALLDEF(phylo_phy_node_rotate, 2),
    CALLDEF(phylo_plot_cartesian, 3),
    CALLDEF(phylo_plot_polar, 2),
    {NULL, NULL, 0}
};


void attribute_visible R_init_phylo(DllInfo *info)
{
    R_registerRoutines(info, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(info, FALSE);
    R_forceSymbols(info, TRUE);

    /* Register the C APIs for use by other packages */
    R_RegisterCCallable(
        "phylo", "phy_node_alloc", (DL_FUNC) &phy_node_alloc);
    R_RegisterCCallable(
        "phylo", "phy_node_free", (DL_FUNC) &phy_node_free);
    R_RegisterCCallable(
        "phylo", "phy_node_add_child", (DL_FUNC) &phy_node_add_child);
    R_RegisterCCallable(
        "phylo", "phy_node_prune", (DL_FUNC) &phy_node_prune);
    R_RegisterCCallable(
        "phylo", "phy_build", (DL_FUNC) &phy_build);
    R_RegisterCCallable(
        "phylo", "phy_read_newickstr", (DL_FUNC) &phy_read_newickstr);
    R_RegisterCCallable(
        "phylo", "phy_write_newickstr", (DL_FUNC) &phy_write_newickstr);
    R_RegisterCCallable(
        "phylo", "phy_read_newickfile", (DL_FUNC) &phy_read_newickfile);
    R_RegisterCCallable(
        "phylo", "phy_write_newickfile", (DL_FUNC) &phy_write_newickfile);
    R_RegisterCCallable(
        "phylo", "phy_free", (DL_FUNC) &phy_free);
    R_RegisterCCallable(
        "phylo", "phy_cursor_prepare", (DL_FUNC) &phy_cursor_prepare);
    R_RegisterCCallable(
        "phylo", "phy_cursor_prepare_v2", (DL_FUNC) &phy_cursor_prepare_v2);
    R_RegisterCCallable(
        "phylo", "phy_cursor_step", (DL_FUNC) &phy_cursor_step);
    R_RegisterCCallable(
        "phylo", "phy_isbinary", (DL_FUNC) &phy_isbinary);
    R_RegisterCCallable(
        "phylo", "phy_extract_clade", (DL_FUNC) &phy_extract_clade);
    R_RegisterCCallable(
        "phylo", "phy_extract_subtree", (DL_FUNC) &phy_extract_subtree);
    R_RegisterCCallable(
        "phylo", "phy_node_rotate", (DL_FUNC) &phy_node_rotate);
    R_RegisterCCallable(
        "phylo", "phy_duplicate", (DL_FUNC) &phy_duplicate);
    R_RegisterCCallable(
        "phylo", "phy_reroot", (DL_FUNC) &phy_reroot);
    R_RegisterCCallable(
        "phylo", "phy_unroot", (DL_FUNC) &phy_unroot);
    R_RegisterCCallable(
        "phylo", "phy_node_brlen", (DL_FUNC) &phy_node_brlen);
    R_RegisterCCallable(
        "phylo", "phy_node_index", (DL_FUNC) &phy_node_index);
    R_RegisterCCallable(
        "phylo", "phy_node_ndesc", (DL_FUNC) &phy_node_ndesc);
    R_RegisterCCallable(
        "phylo", "phy_node_istip", (DL_FUNC) &phy_node_istip);
    R_RegisterCCallable(
        "phylo", "phy_node_lfdesc", (DL_FUNC) &phy_node_lfdesc);
    R_RegisterCCallable(
        "phylo", "phy_node_rtdesc", (DL_FUNC) &phy_node_rtdesc);
    R_RegisterCCallable(
        "phylo", "phy_node_anc", (DL_FUNC) &phy_node_anc);
    R_RegisterCCallable(
        "phylo", "phy_node_next", (DL_FUNC) &phy_node_next);
    R_RegisterCCallable(
        "phylo", "phy_node_prev", (DL_FUNC) &phy_node_prev);
    R_RegisterCCallable(
        "phylo", "phy_isrooted", (DL_FUNC) &phy_isrooted);
    R_RegisterCCallable(
        "phylo", "phy_nnode", (DL_FUNC) &phy_nnode);
    R_RegisterCCallable(
        "phylo", "phy_ntip", (DL_FUNC) &phy_ntip);
    R_RegisterCCallable(
        "phylo", "phy_node_swap", (DL_FUNC) &phy_node_swap);
    R_RegisterCCallable(
        "phylo", "phy_ladderize", (DL_FUNC) &phy_ladderize);
    R_RegisterCCallable(
        "phylo", "phy_root", (DL_FUNC) &phy_root);
    R_RegisterCCallable(
        "phylo", "phy_node_get", (DL_FUNC) &phy_node_get);
    R_RegisterCCallable(
        "phylo", "phy_node_find", (DL_FUNC) &phy_node_find);
    R_RegisterCCallable(
        "phylo", "phy_node_set_data", (DL_FUNC) &phy_node_set_data);
    R_RegisterCCallable(
        "phylo", "phy_node_set_index", (DL_FUNC) &phy_node_set_index);
    R_RegisterCCallable(
        "phylo", "phy_node_set_brlen", (DL_FUNC) &phy_node_set_brlen);
    R_RegisterCCallable(
        "phylo", "phy_node_set_label", (DL_FUNC) &phy_node_set_label);
    R_RegisterCCallable(
        "phylo", "phy_node_data", (DL_FUNC) &phy_node_data);
    R_RegisterCCallable(
        "phylo", "phy_node_spanning_pair", (DL_FUNC) &phy_node_spanning_pair);
    R_RegisterCCallable(
        "phylo", "phy_node_spanning_index", (DL_FUNC) &phy_node_spanning_index);
    R_RegisterCCallable(
        "phylo", "phy_node_mrca", (DL_FUNC) &phy_node_mrca);
    R_RegisterCCallable(
        "phylo", "phy_node_foreach", (DL_FUNC) &phy_node_foreach);
    R_RegisterCCallable(
        "phylo", "phy_cursor_alloc", (DL_FUNC) &phy_cursor_alloc);
    R_RegisterCCallable(
        "phylo", "phy_cursor_free", (DL_FUNC) &phy_cursor_free);
    R_RegisterCCallable(
        "phylo", "phy_node_label", (DL_FUNC) &phy_node_label);
    R_RegisterCCallable(
        "phylo", "phy_node_note", (DL_FUNC) &phy_node_note);
    R_RegisterCCallable(
        "phylo", "phy_errmsg", (DL_FUNC) &phy_errmsg);
}
