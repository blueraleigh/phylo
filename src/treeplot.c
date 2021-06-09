#include <R.h>
#include <Rinternals.h>
#include "phy.h"


SEXP phylo_plot_cartesian(SEXP rtree, SEXP ages, SEXP direction)
{
    struct phy *phy = (struct phy *)R_ExternalPtrAddr(rtree);
    int i;
    int d = INTEGER(direction)[0];
    int nnode = phy_nnode(phy);
    int ntip = phy_ntip(phy);

    SEXP coord = PROTECT(allocMatrix(REALSXP, nnode, 4));
    SEXP bar = PROTECT(allocMatrix(REALSXP, (nnode-ntip),  4));

    double a;
    double b;
    double *segs = REAL(coord);
    double *bars = REAL(bar);
    double *age = REAL(ages);
    struct phy_node *lf;
    struct phy_node *rt;
    struct phy_node *node;
    struct phy_cursor *cursor;

    double ypos = ntip;
    double maxage = 0;

    for (i = 0; i < nnode; ++i)
    {
        if (age[i] > maxage)
            maxage = age[i];
    }

    cursor = phy_cursor_prepare(phy, phy_root(phy), ALL_NODES, POSTORDER);

    while ((node = phy_cursor_step(cursor)) != 0)
    {
        i = phy_node_index(node);
        lf = phy_node_lfdesc(node);
        rt = phy_node_rtdesc(node);
        switch (d)
        {
            case 0:
                segs[i + 0 * nnode] = age[i];                        // x0
                segs[i + 1 * nnode] = age[i] - phy_node_brlen(node); // x1
                if (phy_node_istip(node))
                {
                    segs[i + 2 * nnode] = ypos;                 // y0
                    segs[i + 3 * nnode] = ypos--;               // y1
                }
                else
                {
                    a = segs[phy_node_index(lf) + 2 * nnode];
                    b = segs[phy_node_index(rt) + 2 * nnode];
                    segs[i + 2 * nnode] = (a + b) / 2;
                    segs[i + 3 * nnode] = (a + b) / 2;
                    bars[(i-ntip) + 0 * (nnode-ntip)] = age[i];
                    bars[(i-ntip) + 1 * (nnode-ntip)] = age[i];
                    bars[(i-ntip) + 2 * (nnode-ntip)] = a;
                    bars[(i-ntip) + 3 * (nnode-ntip)] = b;
                }
                break;
            case 1:
                segs[i + 0 * nnode] = maxage - age[i];
                segs[i + 1 * nnode] = maxage - age[i] + phy_node_brlen(node);
                if (phy_node_istip(node))
                {
                    segs[i + 2 * nnode] = ypos;
                    segs[i + 3 * nnode] = ypos--;
                }
                else
                {
                    a = segs[phy_node_index(lf) + 2 * nnode];
                    b = segs[phy_node_index(rt) + 2 * nnode];
                    segs[i + 2 * nnode] = (a + b) / 2;
                    segs[i + 3 * nnode] = (a + b) / 2;
                    bars[(i-ntip) + 0 * (nnode-ntip)] = maxage - age[i];
                    bars[(i-ntip) + 1 * (nnode-ntip)] = maxage - age[i];
                    bars[(i-ntip) + 2 * (nnode-ntip)] = a;
                    bars[(i-ntip) + 3 * (nnode-ntip)] = b;
                }
                break;
            case 2:
                segs[i + 2 * nnode] = age[i];
                segs[i + 3 * nnode] = age[i] - phy_node_brlen(node);
                if (phy_node_istip(node))
                {
                    segs[i + 0 * nnode] = ypos;
                    segs[i + 1 * nnode] = ypos--;
                }
                else
                {
                    a = segs[phy_node_index(lf) + 0 * nnode];
                    b = segs[phy_node_index(rt) + 0 * nnode];
                    segs[i + 0 * nnode] = (a + b) / 2;
                    segs[i + 1 * nnode] = (a + b) / 2;
                    bars[(i-ntip) + 2 * (nnode-ntip)] = age[i];
                    bars[(i-ntip) + 3 * (nnode-ntip)] = age[i];
                    bars[(i-ntip) + 0 * (nnode-ntip)] = a;
                    bars[(i-ntip) + 1 * (nnode-ntip)] = b;
                }
                break;
            case 3:
                segs[i + 2 * nnode] = maxage - age[i];
                segs[i + 3 * nnode] = maxage - age[i] + phy_node_brlen(node);
                if (phy_node_istip(node))
                {
                    segs[i + 0 * nnode] = ypos;
                    segs[i + 1 * nnode] = ypos--;
                }
                else
                {
                    a = segs[phy_node_index(lf) + 0 * nnode];
                    b = segs[phy_node_index(rt) + 0 * nnode];
                    segs[i + 0 * nnode] = (a + b) / 2;
                    segs[i + 1 * nnode] = (a + b) / 2;
                    bars[(i-ntip) + 2 * (nnode-ntip)] = maxage - age[i];
                    bars[(i-ntip) + 3 * (nnode-ntip)] = maxage - age[i];
                    bars[(i-ntip) + 0 * (nnode-ntip)] = a;
                    bars[(i-ntip) + 1 * (nnode-ntip)] = b;
                }
                break;
        }
    }
    SEXP ret = PROTECT(allocVector(VECSXP, 2));
    SET_VECTOR_ELT(ret, 0, coord);
    SET_VECTOR_ELT(ret, 1, bar);
    UNPROTECT(3);
    return ret;
}


SEXP phylo_plot_polar(SEXP rtree, SEXP step)
{
    struct phy *phy = (struct phy *)R_ExternalPtrAddr(rtree);
    int i;
    int z = 0;
    int nnode = phy_nnode(phy);
    SEXP theta = PROTECT(allocMatrix(REALSXP, nnode, 3));
    double *th = REAL(theta);
    double a, b, vstep = REAL(step)[0];
    struct phy_node *lf;
    struct phy_node *rt;
    struct phy_node *node;
    struct phy_cursor *cursor;

    cursor = phy_cursor_prepare(phy, phy_root(phy), ALL_NODES, POSTORDER);

    while ((node = phy_cursor_step(cursor)) != 0)
    {
        i = phy_node_index(node);
        lf = phy_node_lfdesc(node);
        rt = phy_node_rtdesc(node);
        if (phy_node_istip(node))
        {
            th[i + 0 * nnode] = vstep * z++;
            th[i + 1 * nnode] = 0;
            th[i + 2 * nnode] = 0;
        }
        else
        {
            a = th[phy_node_index(lf)];
            b = th[phy_node_index(rt)];
            th[i + 0 * nnode] = (a + b) / 2;
            th[i + 1 * nnode] = a;
            th[i + 2 * nnode] = b;
        }
    }
    UNPROTECT(1);
    return theta;
}

