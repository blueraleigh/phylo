/* This file needs to stay in sync with ../../src/phy.h
**
** This file just provides wrappers to all the functions
** in ../../src/phy.h so that they can be called from C
** code in other R packages. If another R package wants
** to make use of these APIs it simply needs to add the
** phylo package as an entry in the Depends and LinkingTo
** targets in its Description file. It can then #include
** this file in its own sources. Be sure to #define
** PHY_API_IMPLEMENTATION before #include'ing (and only in
** one place if there will be more than one #include).
*/

#ifndef PHY_API_H
#define PHY_API_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <R_ext/Rdynload.h>

#ifdef __cplusplus
extern "C" {
#endif

#define LIBPHY_VERSION "1.0.0"

/* Macros used for tree traversals. */
#define PREORDER 0
#define POSTORDER 1
#define ALL_NODES 0
#define INTERNAL_NODES_ONLY 1

#define PHY_OK 0
#define PHY_ERR 1

struct phy;
struct phy_node;
struct phy_cursor;

// Allocate a new node. Return 0 on success, 1 on failure.
int phy_node_alloc(struct phy_node **node);

// Free an allocated node.
void phy_node_free(struct phy_node *node);

// Add q to p's list of immediate descendants.
void phy_node_add_child(struct phy_node *p, struct phy_node *q);

// Remove q from p's list of immediate descendants and return q.
// If q is not a child of p return NULL.
struct phy_node *phy_node_prune(struct phy_node *p, struct phy_node *q);

// Build a phylogeny from the connected nodes rooted at *root.
struct phy *phy_build(struct phy_node *root, int nnode, int ntip);

// Build a phylogeny from a newick string. The returned phy object must be
// free'd with phy_free.
struct phy *phy_read_newickstr(const char *newick);

// Write a phylogeny to a newick string
char *phy_write_newickstr(struct phy *phy);

// Build a phylogeny from a newick file. The returned phy object must be
// free'd with phy_free.
struct phy *phy_read_newickfile(const char *filename);

// Write a phylogeny to a newick file. Returns 1 on error, 0 on success
int phy_write_newickfile(
    struct phy *phy, const char *filename, const char *mode);

// Free memory allocated to a phylogeny
void phy_free(struct phy *phy);

/* The canonical way to perform a tree traversal from
** a particular node is like so,
**
** struct phy_cursor *cursor = 0;
** cursor = phy_cursor_prepare(phy, node, ALL_NODES, PREORDER);
** while ((node = phy_cursor_step(cursor)) != 0)
** {
**      // perform some operations on node
** }
**
** Which will visit all nodes in the tree in preorder
** traversal sequence. Other permissible values are,
**
**      phy_cursor_prepare(phy, node, ALL_NODES, POSTORDER);
**      phy_cursor_prepare(phy, node, INTERNAL_NODES_ONLY, PREORDER);
**      phy_cursor_prepare(phy, node, INTERNAL_NODES_ONLY, POSTORDER);
**
** Attempting to perform a traversal by using the
** loop construct without calling phy_cursor_prepare
** will break from the loop immediately. */

// Prepare a cursor for phylogeny traversal. In this implementation the
// cursor that is returned does not need to be freed. This version is
// not thread safe. All calls to this function return the same cursor
// object. If multiple threads need to be traversing the phylogeny
// simultaneously they should use phy_cursor_prepare_v2 in conjunction
// with independent cursor objects obtained via phy_cursor_alloc.
struct phy_cursor *phy_cursor_prepare(
    struct phy *phy,
    struct phy_node *node,
    int visit,
    int order
);

// Users should arrange to free the cursor with phy_cursor_free when it is
// no longer needed.
void phy_cursor_prepare_v2(
    struct phy *phy,
    struct phy_node *node,
    struct phy_cursor *cursor,
    int visit,
    int order
);

// Advance the cursor to the next node
struct phy_node *phy_cursor_step(struct phy_cursor *cursor);

// Test whether a phylogeny is bifurcating
int phy_isbinary(struct phy *phy);

// Return the subtree rooted at a given node.
struct phy *phy_extract_clade(struct phy_node *node);

// Return the connected subtree defined by a set of terminal nodes
struct phy *phy_extract_subtree(
    int ntip, struct phy_node **tips, struct phy *phy);

// Rotate a set of nodes
void phy_node_rotate(int n, struct phy_node **nodes, struct phy *phy);

// Return a deep copy of a phylogeny.
struct phy *phy_duplicate(struct phy *phy);

// Re-root the *in phylogeny on node, storing the re-rooted tree
// in *out. If *in is equal to *out, *in is first free'd.
void phy_reroot(struct phy_node *node, struct phy **in, struct phy **out);

// Un-root the root phylogeny *in and store it in *out. If *in
// is already unrooted sets *out equal to NULL.
// If *in is equal to *out, *in is first free'd.
void phy_unroot(struct phy **in, struct phy **out);

// Return the branch length subtending a node
double phy_node_brlen(struct phy_node *node);

// Return the index of a node
int phy_node_index(struct phy_node *node);

// Return the number of immediate descendants of a node
int phy_node_ndesc(struct phy_node *node);

// Test whether a node is terminal
int phy_node_istip(struct phy_node *node);

// Return a node's immediate leftmost descendant
struct phy_node *phy_node_lfdesc(struct phy_node *node);

// Return a node's immediate rightmost descendant
struct phy_node *phy_node_rtdesc(struct phy_node *node);

// Return a node's immediate ancestor
struct phy_node *phy_node_anc(struct phy_node *node);

// Return a node's next sibling
struct phy_node *phy_node_next(struct phy_node *node);

// Return a node's previous sibling
struct phy_node *phy_node_prev(struct phy_node *node);

/* For non-bifurcating trees, to loop over all the
** immediate descendants of a node use
**
**   struct phy_node *d;
**   for (d = phy_node_lfdesc(node); d != 0; d = phy_node_next(d))
**   {
**      // do something with descendant
**   }
**
*/

// Determine if a phylogeny is rooted. If there is a basal polytomy the
// tree is considered unrooted and 0 is returned. Otherwise, 1.
int phy_isrooted(struct phy *phy);

// Return the number of nodes in a phylogeny
int phy_nnode(struct phy *phy);

// Return the number of terminal nodes in a phylogeny
int phy_ntip(struct phy *phy);

// Swap the position of a and b in the child list
void phy_node_swap(struct phy_node *a, struct phy_node *b);

// Ladderize the phylogeny
void phy_ladderize(struct phy *phy, int *n, int *perm);

// Return the root node of a phylogeny
struct phy_node *phy_root(struct phy *phy);

// Return the node with the given index from a phylogeny (or NULL if the
// index is invalid).
struct phy_node *phy_node_get(struct phy *phy, int index);

// Return the node with the given label from a phylogeny (or NULL if the
// label is not found).
struct phy_node *phy_node_find(struct phy *phy, const char *label);

// Attach arbitrary client data to a node, passing a destructor (may be NULL).
// Any previously added data is removed (and free'd).
void phy_node_set_data(
    struct phy_node *node, void *data, void (*data_free)(void *));

// Set the index for a node.
void phy_node_set_index(struct phy_node *node, int index);

// Set the length of the branch subtending a node.
void phy_node_set_brlen(
    struct phy_node *node, double brlen);

// Set the label for a node.
void phy_node_set_label(
    struct phy_node *node, const char *label);

// Return client data attached to a node
void *phy_node_data(struct phy_node *node);

// Return the labels of the most distant terminal node's whose most
// recent common ancestor is *node and store them in *a and *b.
void phy_node_spanning_pair(
    struct phy_node *node, const char **a, const char **b);

// Return the indices of the most distant terminal node's whose most
// recent common ancestor is *node and store them in *a and *b.
void phy_node_spanning_index(
    struct phy_node *node, int *a, int *b);

// Return the most recent common ancestor of a and b.
struct phy_node *phy_node_mrca(
    struct phy *phy, struct phy_node *a, struct phy_node *b);

// Apply function FUN to each node visited by a specified type of tree
// traversal.
void phy_node_foreach(
    struct phy *phy,
    struct phy_node *node,
    int visit,
    int order,
    void (*FUN)(struct phy_node *node, struct phy *phy, void *param),
    void *param);

// Allocate memory for a new phy_cursor object. Return 0 on success; 1, error
int phy_cursor_alloc(struct phy_cursor **cursor);

// Free memory allocated to a phy_cursor object
void phy_cursor_free(struct phy_cursor *cursor);

const char *phy_node_label(struct phy_node *node);

const char *phy_node_note(struct phy_node *node);

// Return the current error message
const char *phy_errmsg();


#ifdef PHY_API_IMPLEMENTATION

int phy_node_alloc(struct phy_node **node)
{
    static int(*fun)(struct phy_node **) = NULL;
    if (!fun)
    {
        fun = (int(*)(struct phy_node **))R_GetCCallable(
            "phylo", "phy_node_alloc");
    }
    return fun(node);
}

void phy_node_free(struct phy_node *node)
{
    static void(*fun)(struct phy_node *node) = NULL;
    if (!fun)
    {
        fun = (void(*)(struct phy_node *))R_GetCCallable(
            "phylo", "phy_node_free");
    }
    fun(node);
}

void phy_node_add_child(struct phy_node *p, struct phy_node *q)
{
    static void(*fun)(struct phy_node *, struct phy_node *) = NULL;
    if (!fun)
    {
        fun = (void(*)(struct phy_node *, struct phy_node *))R_GetCCallable(
            "phylo", "phy_node_add_child");
    }
    fun(p, q);
}

struct phy_node *phy_node_prune(struct phy_node *p, struct phy_node *q)
{
    static struct phy_node *(*fun)(struct phy_node *, struct phy_node *) = NULL;
    if (!fun)
    {
        fun = (struct phy_node *(*)(struct phy_node *, struct phy_node *))
            R_GetCCallable("phylo", "phy_node_prune");
    }
    return fun(p, q);
}

struct phy *phy_build(struct phy_node *root, int nnode, int ntip)
{
    static struct phy *(*fun)(struct phy_node *, int, int) = NULL;
    if (!fun)
    {
        fun = (struct phy *(*)(struct phy_node *, int, int))R_GetCCallable(
            "phylo", "phy_build");
    }
    return fun(root, nnode, ntip);
}

struct phy *phy_read_newickstr(const char *newick)
{
    static struct phy *(*fun)(const char *) = NULL;
    if (!fun)
    {
        fun = (struct phy *(*)(const char *))R_GetCCallable(
            "phylo", "phy_read_newickstr");
    }
    return fun(newick);
}

char *phy_write_newickstr(struct phy *phy)
{
    static char *(*fun)(struct phy *) = NULL;
    if (!fun)
    {
        fun = (char *(*)(struct phy *))R_GetCCallable(
            "phylo", "phy_write_newickstr");
    }
    return fun(phy);
}

struct phy *phy_read_newickfile(const char *filename)
{
    static struct phy *(*fun)(const char *) = NULL;
    if (!fun)
    {
        fun = (struct phy *(*)(const char *))R_GetCCallable(
            "phylo", "phy_read_newickfile");
    }
    return fun(filename);
}

int phy_write_newickfile(
    struct phy *phy, const char *filename, const char *mode)
{
    static int(*fun)(struct phy *, const char *, const char *) = NULL;
    if (!fun)
    {
        fun = (int(*)(struct phy *, const char *, const char *))R_GetCCallable(
            "phylo", "phy_write_newickfile");
    }
    return fun(phy, filename, mode);
}

void phy_free(struct phy *phy)
{
    static void(*fun)(struct phy *) = NULL;
    if (!fun)
    {
        fun = (void(*)(struct phy *))R_GetCCallable("phylo", "phy_free");
    }
    fun(phy);
}

struct phy_cursor *phy_cursor_prepare(
    struct phy *phy,
    struct phy_node *node,
    int visit,
    int order
)
{
    static struct phy_cursor *(*fun)(
        struct phy *, struct phy_node *, int, int) = NULL;
    if (!fun)
    {
        fun = (struct phy_cursor *(*)(struct phy *, struct phy_node *, int, int))
            R_GetCCallable("phylo", "phy_cursor_prepare");
    }
    return fun(phy, node, visit, order);
}


void phy_cursor_prepare_v2(
    struct phy *phy,
    struct phy_node *node,
    struct phy_cursor *cursor,
    int visit,
    int order
)
{
    static void(*fun)(
        struct phy *, struct phy_node *, struct phy_cursor *, int, int) = NULL;
    if (!fun)
    {
        fun = (void(*)(struct phy *, struct phy_node *, struct phy_cursor *, int, int))
            R_GetCCallable("phylo", "phy_cursor_prepare_v2");
    }
    fun(phy, node, cursor, visit, order);
}

struct phy_node *phy_cursor_step(struct phy_cursor *cursor)
{
    static struct phy_node *(*fun)(struct phy_cursor *) = NULL;
    if (!fun)
    {
        fun = (struct phy_node *(*)(struct phy_cursor *))R_GetCCallable(
            "phylo", "phy_cursor_step");
    }
    return fun(cursor);
}

int phy_isbinary(struct phy *phy)
{
    static int(*fun)(struct phy *) = NULL;
    if (!fun)
    {
        fun = (int(*)(struct phy *))R_GetCCallable("phylo", "phy_isbinary");
    }
    return fun(phy);
}

struct phy *phy_extract_clade(struct phy_node *node)
{
    static struct phy *(*fun)(struct phy_node *) = NULL;
    if (!fun)
    {
        fun = (struct phy *(*)(struct phy_node *))R_GetCCallable(
            "phylo", "phy_extract_clade");
    }
    return fun(node);
}

struct phy *phy_extract_subtree(
    int ntip, struct phy_node **tips, struct phy *phy)
{
    static struct phy *(*fun)(int, struct phy_node **, struct phy *) = NULL;
    if (!fun)
    {
        fun = (struct phy *(*)(int, struct phy_node **, struct phy *))R_GetCCallable(
            "phylo", "phy_extract_subtree");
    }
    return fun(ntip, tips, phy);
}

void phy_node_rotate(int n, struct phy_node **nodes, struct phy *phy)
{
    static void(*fun)(int, struct phy_node **, struct phy *) = NULL;
    if (!fun)
    {
        fun = (void(*)(int, struct phy_node **, struct phy *))R_GetCCallable(
            "phylo", "phy_node_rotate");
    }
    fun(n, nodes, phy);
}

struct phy *phy_duplicate(struct phy *phy)
{
    static struct phy *(*fun)(struct phy *) = NULL;
    if (!fun)
    {
        fun = (struct phy *(*)(struct phy *))R_GetCCallable(
            "phylo", "phy_duplicate");
    }
    return fun(phy);
}

void phy_reroot(struct phy_node *node, struct phy **in, struct phy **out)
{
    static void(*fun)(struct phy_node *, struct phy **, struct phy **) = NULL;
    if (!fun)
    {
        fun = (void(*)(struct phy_node *, struct phy **, struct phy **))
            R_GetCCallable("phylo", "phy_reroot");
    }
    fun(node, in, out);
}

void phy_unroot(struct phy **in, struct phy **out)
{
    static void(*fun)(struct phy **, struct phy **) = NULL;
    if (!fun)
    {
        fun = (void(*)(struct phy **, struct phy **))R_GetCCallable(
            "phylo", "phy_unroot");
    }
    fun(in, out);
}

double phy_node_brlen(struct phy_node *node)
{
    static double(*fun)(struct phy_node *) = NULL;
    if (!fun)
    {
        fun = (double(*)(struct phy_node *))R_GetCCallable(
            "phylo", "phy_node_brlen");
    }
    return fun(node);
}

int phy_node_index(struct phy_node *node)
{
    static int(*fun)(struct phy_node *) = NULL;
    if (!fun)
    {
        fun = (int(*)(struct phy_node *))R_GetCCallable(
            "phylo", "phy_node_index");
    }
    return fun(node);
}

int phy_node_ndesc(struct phy_node *node)
{
    static int(*fun)(struct phy_node *) = NULL;
    if (!fun)
    {
        fun = (int(*)(struct phy_node *))R_GetCCallable(
            "phylo", "phy_node_ndesc");
    }
    return fun(node);
}

int phy_node_istip(struct phy_node *node)
{
    static int(*fun)(struct phy_node *) = NULL;
    if (!fun)
    {
        fun = (int(*)(struct phy_node *))R_GetCCallable(
            "phylo", "phy_node_istip");
    }
    return fun(node);
}

struct phy_node *phy_node_lfdesc(struct phy_node *node)
{
    static struct phy_node *(*fun)(struct phy_node *) = NULL;
    if (!fun)
    {
        fun = (struct phy_node *(*)(struct phy_node *))R_GetCCallable(
            "phylo", "phy_node_lfdesc");
    }
    return fun(node);
}

struct phy_node *phy_node_rtdesc(struct phy_node *node)
{
    static struct phy_node *(*fun)(struct phy_node *) = NULL;
    if (!fun)
    {
        fun = (struct phy_node *(*)(struct phy_node *))R_GetCCallable(
            "phylo", "phy_node_rtdesc");
    }
    return fun(node);
}

struct phy_node *phy_node_anc(struct phy_node *node)
{
    static struct phy_node *(*fun)(struct phy_node *) = NULL;
    if (!fun)
    {
        fun = (struct phy_node *(*)(struct phy_node *))R_GetCCallable(
            "phylo", "phy_node_anc");
    }
    return fun(node);
}

struct phy_node *phy_node_next(struct phy_node *node)
{
    static struct phy_node *(*fun)(struct phy_node *) = NULL;
    if (!fun)
    {
        fun = (struct phy_node *(*)(struct phy_node *))R_GetCCallable(
            "phylo", "phy_node_next");
    }
    return fun(node);
}

struct phy_node *phy_node_prev(struct phy_node *node)
{
    static struct phy_node *(*fun)(struct phy_node *) = NULL;
    if (!fun)
    {
        fun = (struct phy_node *(*)(struct phy_node *))R_GetCCallable(
            "phylo", "phy_node_prev");
    }
    return fun(node);
}

int phy_isrooted(struct phy *phy)
{
    static int(*fun)(struct phy *) = NULL;
    if (!fun)
    {
        fun = (int(*)(struct phy *))R_GetCCallable("phylo", "phy_isrooted");
    }
    return fun(phy);
}

int phy_nnode(struct phy *phy)
{
    static int(*fun)(struct phy *) = NULL;
    if (!fun)
    {
        fun = (int(*)(struct phy *))R_GetCCallable("phylo", "phy_nnode");
    }
    return fun(phy);
}

int phy_ntip(struct phy *phy)
{
    static int(*fun)(struct phy *) = NULL;
    if (!fun)
    {
        fun = (int(*)(struct phy *))R_GetCCallable("phylo", "phy_ntip");
    }
    return fun(phy);
}

void phy_node_swap(struct phy_node *a, struct phy_node *b)
{
    static void(*fun)(struct phy_node *, struct phy_node *) = NULL;
    if (!fun)
    {
        fun = (void(*)(struct phy_node *, struct phy_node *))R_GetCCallable(
            "phylo", "phy_node_swap");
    }
    fun(a, b);
}

void phy_ladderize(struct phy *phy, int *n, int *perm)
{
    static void(*fun)(struct phy *, int *, int *) = NULL;
    if (!fun)
    {
        fun = (void(*)(struct phy *, int *, int *))R_GetCCallable(
            "phylo", "phy_ladderize");
    }
    fun(phy, n, perm);
}

struct phy_node *phy_root(struct phy *phy)
{
    static struct phy_node *(*fun)(struct phy *) = NULL;
    if (!fun)
    {
        fun = (struct phy_node *(*)(struct phy *))R_GetCCallable("phylo", "phy_root");
    }
    return fun(phy);
}

struct phy_node *phy_node_get(struct phy *phy, int index)
{
    static struct phy_node *(*fun)(struct phy *, int) = NULL;
    if (!fun)
    {
        fun = (struct phy_node *(*)(struct phy *, int))R_GetCCallable(
            "phylo", "phy_node_get");
    }
    return fun(phy, index);
}

struct phy_node *phy_node_find(struct phy *phy, const char *label)
{
    static struct phy_node *(*fun)(struct phy *, const char *) = NULL;
    if (!fun)
    {
        fun = (struct phy_node *(*)(struct phy *, const char *))R_GetCCallable(
            "phylo", "phy_node_find");
    }
    return fun(phy, label);
}

void phy_node_set_data(
    struct phy_node *node, void *data, void (*data_free)(void *))
{
    static void (*fun)(struct phy_node *, void *, void (*)(void *)) = NULL;
    if (!fun)
    {
        fun = (void(*)(struct phy_node *, void *, void (*)(void *)))R_GetCCallable(
            "phylo", "phy_node_set_data");
    }
    fun(node, data, data_free);
}

void phy_node_set_index(struct phy_node *node, int index)
{
    static void (*fun)(struct phy_node *, int) = NULL;
    if (!fun)
    {
        fun = (void(*)(struct phy_node *, int))R_GetCCallable(
            "phylo", "phy_node_set_index");
    }
    fun(node, index);
}

void phy_node_set_brlen(struct phy_node *node, double brlen)
{
    static void (*fun)(struct phy_node *, double) = NULL;
    if (!fun)
    {
        fun = (void(*)(struct phy_node *, double))R_GetCCallable(
            "phylo", "phy_node_set_brlen");
    }
    fun(node, brlen);
}

void phy_node_set_label(struct phy_node *node, const char *label)
{
    static void (*fun)(struct phy_node *, const char *) = NULL;
    if (!fun)
    {
        fun = (void(*)(struct phy_node *, const char *))R_GetCCallable(
            "phylo", "phy_node_set_label");
    }
    fun(node, label);
}

void *phy_node_data(struct phy_node *node)
{
    static void *(*fun)(struct phy_node *) = NULL;
    if (!fun)
    {
        fun = (void *(*)(struct phy_node *))R_GetCCallable(
            "phylo", "phy_node_data");
    }
    return fun(node);
}

void phy_node_spanning_pair(
    struct phy_node *node, const char **a, const char **b)
{
    static void(*fun)(struct phy_node *, const char **, const char **) = NULL;
    if (!fun)
    {
        fun = (void(*)(struct phy_node *, const char **, const char **))
            R_GetCCallable("phylo", "phy_node_spanning_pair");
    }
    fun(node, a, b);
}

void phy_node_spanning_index(
    struct phy_node *node, int *a, int *b)
{
    static void(*fun)(struct phy_node *, int *, int *) = NULL;
    if (!fun)
    {
        fun = (void(*)(struct phy_node *, int *, int *))
            R_GetCCallable("phylo", "phy_node_spanning_index");
    }
    fun(node, a, b);
}


struct phy_node *phy_node_mrca(
    struct phy *phy, struct phy_node *a, struct phy_node *b)
{
    static struct phy_node *(*fun)(
        struct phy *, struct phy_node *, struct phy_node *) = NULL;
    if (!fun)
    {
        fun = (struct phy_node *(*)(
            struct phy *, struct phy_node *, struct phy_node *))
        R_GetCCallable("phylo", "phy_node_mrca");
    }
    return fun(phy, a, b);
}

void phy_node_foreach(
    struct phy *phy,
    struct phy_node *node,
    int visit,
    int order,
    void (*FUN)(struct phy_node *node, struct phy *phy, void *param),
    void *param)
{
    static void(*fun)(
        struct phy *,
        struct phy_node *,
        int,
        int,
        void (*)(struct phy_node *, struct phy *, void *),
        void *) = NULL;
    if (!fun)
    {
        fun = (void(*)(
            struct phy *,
            struct phy_node *,
            int,
            int,
            void (*)(struct phy_node *, struct phy *, void *),
            void *))
        R_GetCCallable("phylo", "phy_node_foreach");
    }
    fun(phy, node, visit, order, FUN, param);
}

int phy_cursor_alloc(struct phy_cursor **cursor)
{
    static int(*fun)(struct phy_cursor **) = NULL;
    if (!fun)
    {
        fun = (int(*)(struct phy_cursor **))R_GetCCallable(
                "phylo", "phy_cursor_alloc");
    }
    return fun(cursor);
}

void phy_cursor_free(struct phy_cursor *cursor)
{
    static void(*fun)(struct phy_cursor *) = NULL;
    if (!fun)
    {
        fun = (void(*)(struct phy_cursor *))R_GetCCallable(
            "phylo", "phy_cursor_free");
    }
    fun(cursor);
}

const char *phy_node_label(struct phy_node *node)
{
    static const char *(*fun)(struct phy_node *) = NULL;
    if (!fun)
    {
        fun = (const char *(*)(struct phy_node *))R_GetCCallable(
            "phylo", "phy_node_label");
    }
    return fun(node);
}

const char *phy_node_note(struct phy_node *node)
{
    static const char *(*fun)(struct phy_node *) = NULL;
    if (!fun)
    {
        fun = (const char *(*)(struct phy_node *))R_GetCCallable(
            "phylo", "phy_node_note");
    }
    return fun(node);
}

const char *phy_errmsg()
{
    static const char *(*fun)() = NULL;
    if (!fun)
    {
        fun = (const char *(*)())R_GetCCallable(
            "phylo", "phy_errmsg");
    }
    return fun();
}

#endif /* PHY_API_IMPLEMENTATION */

#ifdef __cplusplus
}
#endif

#endif /* PHY_API_H */
