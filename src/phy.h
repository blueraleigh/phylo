#ifndef PHY_H
#define PHY_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>


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
void phy_ladderize(struct phy *phy, int *perm);

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
void phy_node_add_data(
    struct phy_node *node, void *data, void (*data_free)(void *));

// Return client data attached to a node
void *phy_node_data(struct phy_node *node);

// Return the labels of the most distant terminal node's whose most
// recent common ancestor is *node and store them in *a and *b.
void phy_node_spanning_pair(
    struct phy_node *node, const char **a, const char **b);

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

#ifdef __cplusplus
}
#endif

#endif /* PHY_H */
