#include <R.h>
#include <Rinternals.h>
#include "phy.h"


void phylo_phy_free(SEXP rtree)
{
    if (TYPEOF(rtree) == NILSXP)
        return;
    struct phy *phy = (struct phy*)R_ExternalPtrAddr(rtree);
    phy_free(phy);
    R_ClearExternalPtr(rtree);
}


SEXP phylo_phy_read_newickstr(SEXP newick)
{
    SEXP rtree;
    struct phy *phy = phy_read_newickstr(CHAR(STRING_ELT(newick, 0)));
    if (phy) {
        rtree = PROTECT(R_MakeExternalPtr(phy, R_NilValue, R_NilValue));
        R_RegisterCFinalizer(rtree, &phylo_phy_free);
        setAttrib(rtree, install("root"),
            ScalarInteger(phy_node_index(phy_root(phy))+1));
        setAttrib(rtree, install("Ntip"), ScalarInteger(phy_ntip(phy)));
        setAttrib(rtree, install("Nnode"), ScalarInteger(phy_nnode(phy)));
        UNPROTECT(1);
        return rtree;
    }
    error(phy_errmsg());
}


SEXP phylo_phy_write_newickstr(SEXP rtree)
{
    struct phy *phy = (struct phy*)R_ExternalPtrAddr(rtree);
    char *newick = phy_write_newickstr(phy);
    SEXP ret = PROTECT(allocVector(STRSXP, 1));
    SET_STRING_ELT(ret, 0, mkChar(newick));
    free(newick);
    UNPROTECT(1);
    return ret;
}


SEXP phylo_tiplabels(SEXP rtree)
{
    int i;
    int ntip;
    SEXP tiplabel;
    struct phy *phy = (struct phy*)R_ExternalPtrAddr(rtree);

    ntip = phy_ntip(phy);
    tiplabel = PROTECT(allocVector(STRSXP, ntip));

    for (i = 0; i < ntip; ++i)
    {
        SET_STRING_ELT(tiplabel, i,
            mkChar(phy_node_label(phy_node_get(phy, i))));
    }

    UNPROTECT(1);
    return tiplabel;
}


SEXP phylo_node_notes(SEXP rtree)
{
    int i;
    int nnode;
    const char *note;
    SEXP notes;
    struct phy *phy = (struct phy*)R_ExternalPtrAddr(rtree);

    nnode = phy_nnode(phy);
    notes = PROTECT(allocVector(STRSXP, nnode));

    for (i = 0; i < nnode; ++i) {
        note = phy_node_note(phy_node_get(phy, i));
        if (note)
            SET_STRING_ELT(notes, i, mkChar(note));
        else
            SET_STRING_ELT(notes, i, mkChar(""));
    }

    UNPROTECT(1);
    return notes;
}


SEXP phylo_phy_node_brlens(SEXP rtree)
{
    int i;
    int nnode;
    SEXP brlen;
    struct phy *phy = (struct phy*)R_ExternalPtrAddr(rtree);

    nnode = phy_nnode(phy);
    brlen = PROTECT(allocVector(REALSXP, nnode));

    for (i = 0; i < nnode; ++i)
        REAL(brlen)[i] = phy_node_brlen(phy_node_get(phy, i));

    UNPROTECT(1);
    return brlen;
}


SEXP phylo_phy_node_ages(SEXP rtree)
{
    int i;
    int nnode;
    double *node_age;
    SEXP age;
    struct phy *phy = (struct phy*)R_ExternalPtrAddr(rtree);
    struct phy_node *node;

    nnode = phy_nnode(phy);

    age = PROTECT(allocVector(REALSXP, nnode));
    node_age = REAL(age);

    for (i = 0; i < nnode; ++i)
    {
        node = phy_node_get(phy, i);
        node_age[i] = 0.0;
        while (node != NULL) {
            node_age[i] += phy_node_brlen(node);
            node = phy_node_anc(node);
        }
    }

    UNPROTECT(1);
    return age;
}


SEXP phylo_phy_node_ancestors(SEXP rtree, SEXP node)
{
    int i;
    int sz;
    int end = 0;
    int cnt = 0;
    int nanc = 0;
    struct phy *phy = (struct phy*)R_ExternalPtrAddr(rtree);
    struct phy_node *p = phy_node_get(phy, INTEGER(node)[0]-1);

    SEXP buf;
    SEXP root = PROTECT(list1(buf = allocVector(VECSXP, 1000)));
    SEXP tail = root;

    while (p != NULL) {
        nanc++;
        SET_VECTOR_ELT(buf, cnt++, ScalarInteger(phy_node_index(p)+1));
        if (cnt == 1000) {
            tail = SETCDR(tail, list1(buf = allocVector(VECSXP, 1000)));
            cnt = 0;
        }
        p = phy_node_anc(p);
    }

    SEXP ret = PROTECT(allocVector(INTSXP, nanc));

    while (root != R_NilValue) {
        sz = CDR(root) == R_NilValue ? cnt : 1000;
        for (i = 0; i < sz; ++i)
            INTEGER(ret)[end++] = INTEGER(VECTOR_ELT(CAR(root), i))[0];
        root = CDR(root);
    }

    UNPROTECT(2);
    return ret;
}


SEXP phylo_phy_node_children(SEXP rtree, SEXP node)
{
    int i = 0;
    struct phy *phy = (struct phy*)R_ExternalPtrAddr(rtree);
    struct phy_node *p = phy_node_get(phy, INTEGER(node)[0]-1);
    SEXP ret = PROTECT(allocVector(INTSXP, phy_node_ndesc(p)));
    for (p = phy_node_lfdesc(p); p != 0; p = phy_node_next(p))
        INTEGER(ret)[i++] = phy_node_index(p) + 1;
    UNPROTECT(1);
    return ret;
}


SEXP phylo_phy_node_descendants(SEXP rtree, SEXP node, SEXP visit, SEXP order)
{
    int i;
    int sz;
    int end = 0;
    int cnt = 0;
    int ndesc = 0;
    struct phy_node *d;
    struct phy_cursor *cursor;
    struct phy *phy = (struct phy*)R_ExternalPtrAddr(rtree);

    SEXP buf;
    SEXP root = PROTECT(list1(buf = allocVector(VECSXP, 1000)));
    SEXP tail = root;

    cursor = phy_cursor_prepare(phy, phy_node_get(phy, INTEGER(node)[0]-1),
        INTEGER(visit)[0], INTEGER(order)[0]);

    while ((d = phy_cursor_step(cursor)) != 0) {
        ndesc++;
        SET_VECTOR_ELT(buf, cnt++, ScalarInteger(phy_node_index(d)+1));
        if (cnt == 1000) {
            tail = SETCDR(tail, list1(buf = allocVector(VECSXP, 1000)));
            cnt = 0;
        }
    }

    SEXP descendants = PROTECT(allocVector(INTSXP, ndesc));

    while (root != R_NilValue) {
        sz = CDR(root) == R_NilValue ? cnt : 1000;
        for (i = 0; i < sz; ++i)
            INTEGER(descendants)[end++] = INTEGER(VECTOR_ELT(CAR(root), i))[0];
        root = CDR(root);
    }

    UNPROTECT(2);
    return descendants;
}


SEXP phylo_phy_extract_clade(SEXP rtree, SEXP node)
{
    SEXP rclade;
    struct phy *phy = (struct phy *)R_ExternalPtrAddr(rtree);
    struct phy *clade = phy_extract_clade(
        phy_node_get(phy, INTEGER(node)[0]-1));
    if (clade) {
        rclade = PROTECT(R_MakeExternalPtr(clade, R_NilValue, R_NilValue));
        R_RegisterCFinalizer(rclade, &phylo_phy_free);
        setAttrib(rclade, install("root"),
            ScalarInteger(phy_node_index(phy_root(clade))+1));
        setAttrib(rclade, install("Ntip"), ScalarInteger(phy_ntip(clade)));
        setAttrib(rclade, install("Nnode"), ScalarInteger(phy_nnode(clade)));
        UNPROTECT(1);
        return rclade;
    }
    error(phy_errmsg());
}


SEXP phylo_phy_extract_subtree(SEXP rtree, SEXP ntip, SEXP tips)
{
    SEXP rsubtree;
    int i;
    int Ntip = INTEGER(ntip)[0];
    struct phy *phy = (struct phy *)R_ExternalPtrAddr(rtree);
    struct phy_node *nodes[Ntip];

    for (i = 0; i < Ntip; ++i)
        nodes[i] = phy_node_get(phy, INTEGER(tips)[i]-1);

    struct phy *subtree = phy_extract_subtree(Ntip, nodes, phy);
    if (subtree) {
        rsubtree = PROTECT(R_MakeExternalPtr(subtree, R_NilValue, R_NilValue));
        R_RegisterCFinalizer(rsubtree, &phylo_phy_free);
        setAttrib(rsubtree, install("root"),
            ScalarInteger(phy_node_index(phy_root(subtree))+1));
        setAttrib(rsubtree, install("Ntip"), ScalarInteger(phy_ntip(subtree)));
        setAttrib(rsubtree, install("Nnode"), ScalarInteger(phy_nnode(subtree)));
        UNPROTECT(1);
        return rsubtree;
    }
    error(phy_errmsg());
}


SEXP phylo_phy_ladderize(SEXP rtree, SEXP ndesc)
{
    struct phy *phy = (struct phy *)R_ExternalPtrAddr(rtree);

    SEXP perm = PROTECT(allocVector(INTSXP, phy_nnode(phy)));

    phy_ladderize(phy, INTEGER(ndesc), INTEGER(perm));

    for (int i = 0; i < phy_nnode(phy); ++i)
        INTEGER(perm)[i] += 1;

    UNPROTECT(1);
    return perm;
}


SEXP phylo_phy_node_rotate(SEXP rtree, SEXP index)
{
    int i;
    int n = LENGTH(index);
    struct phy *phy;
    struct phy_node *nodes[n];

    phy = (struct phy *)R_ExternalPtrAddr(rtree);


    for (i = 0; i < n; ++i)
        nodes[i] = phy_node_get(phy, INTEGER(index)[i]-1);

    phy_node_rotate(n, nodes, phy);

    return R_NilValue;
}
