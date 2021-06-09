#include "phy.h"

#define PHY_ERR1 "cannot allocate memory"
#define PHY_ERR2 "encountered unexpected character in Newick string node label/branch length"
#define PHY_ERR3 "detected unifurcation in Newick string"
#define PHY_ERR4 "malformed Newick string"

static int phy_errno = 0;

/**********************************************************************
**
** Internal functions and definitions of opaque structures
**
**********************************************************************/

// Node structure for a polytomous tree data model
struct phy_node {
    /* Terminal nodes are numbered 0 to N-1, where
    ** N is the number of terminal nodes in the
    ** phylogeny. Internal nodes are numbered N
    ** up to one less than the number of nodes in
    ** the phylogeny. The root node is always N, and
    ** indices are assigned to nodes in a preorder
    ** traversal sequence. */
    int index;

    /* Number of immediate descendants */
    int ndesc;

    /* Name of node */
    char *lab;

    /* Any note that was passed in from a Newick string */
    char *note;

    /* First immediate descendant. The other descendants
    ** are stored as a doubly-linked list with this descendant
    ** as the head. The next descendant is found as lfdesc->next,
    ** the next as lfdesc->next->next, and so on. */
    struct phy_node *lfdesc;

    /* Next sibling node. Sibling nodes share the same ancestor. */
    struct phy_node *next;

    /* Previous sibling node */
    struct phy_node *prev;

    /* Ancestral node */
    struct phy_node *anc;

    /* This is the last node visited in a preorder traversal of
    ** the subtree rooted at this node (null if this node
    ** is a terminal node). It is always a terminal node if not null. */
    struct phy_node *lastvisit;

    /* The length of the branch that leads to this node */
    double brlen;

    /* Arbitrary data for client programs using this library */
    void *data;

    /* Function pointer to free data held by node */
    void (*data_free)(void *);
};


struct phy_cursor {
    int visit;
    int order;
    int cursor;
    struct phy *phy;
    struct phy_node *begin;
    struct phy_node *end;
    struct phy_node *next;
};


struct phy {
    /* Number of terminal nodes in phylogeny */
    int ntip;

    /* Number of nodes in phylogeny */
    int nnode;

    /* Root node of phylogeny */
    struct phy_node *root;

    /* Array of all nodes arranged in a preorder traversal */
    struct phy_node **nodes;

    /* Array of internal nodes arranged in a preorder traversal.
    ** Note that as a consequence of how nodes are numbered, the
    ** index of an internal node minus the number of tips corresponds
    ** to its position in this array. */
    struct phy_node **inodes;

    /* Preorder visitation sequence of node indices.
    ** For example, vseq[i] will return the position of
    ** the node having index i in the nodes array. */
    int *vseq;

    /* State information for phylogeny traversal */
    struct phy_cursor cursor;
};


struct newick_reader {
    unsigned int n;
    unsigned int nAlloc;
    unsigned int cursor;
    unsigned int ntip;
    unsigned int nnode;
    char *z;
    const char *newick;
    struct phy_node *q;
    struct phy_node *p;
    struct phy_node *root;
};


struct newick_writer {
    unsigned int n;
    unsigned int nAlloc;
    char *newick;
};


/* This function is called on the root node whenever
** an error is encountered during the process of
** building the phylogeny. */
static void cleanup(struct phy_node *root)
{
    struct phy_node *q, *r, *p = root;

    while (p)
    {
        if (p->lfdesc)
            p = p->lfdesc;
        else if (p->next)
            p = p->next;
        else
        {
            // on entry p is a terminal node that marks a clade
            // boundary. this will be the last node visited in a
            // preorder traversal of the subtree rooted at each node
            // on the path back from this node up to and including the
            // first encountered node with a ->next sibling
            q = p->prev;
            while (q != 0)
            {
                r = q;
                q = q->prev;
                free(r->lab);
                free(r->note);
                free(r);
                r = 0;
            }
            while (p->anc != 0 && p->next == 0)
            {
                q = p;
                p = p->anc;
                free(q->lab);
                free(r->note);
                free(q);
                q = 0;
            }
            q = p;
            p = p->next;
            free(q->lab);
            free(r->note);
            free(q);
            q = 0;
        }
    }
}


static struct phy_node *node_new()
{
    struct phy_node *node = malloc(sizeof(struct phy_node));
    if (!node)
    {
        phy_errno = 1;
        return NULL;
    }
    node->index = -1;
    node->ndesc = 0;
    node->lab = 0;
    node->note = 0;
    node->lfdesc = 0;
    node->next = 0;
    node->prev = 0;
    node->anc = 0;
    node->lastvisit = 0;
    node->brlen = 0;
    node->data = 0;
    node->data_free = 0;
    return node;
}


static void node_free(struct phy_node *node)
{
    if (node)
    {
        if (node->data && node->data_free)
            node->data_free(node->data);
        free(node->lab);
        free(node->note);
        free(node);
    }
}


// Read off a label in a Newick string. Return 0 if successful, 1 on error
static int read_label(struct newick_reader *ctx)
{
    char c;
    int toread = 1;
    while (toread)
    {
        c = ctx->newick[ctx->cursor++];
        switch (c)
        {
            case ':':
            case ',':
            case ')':
            case ';':
            case '[':
                ctx->cursor--;
                toread = 0;
                break;
            case ' ':
            case '\n':
            case '\r':
            case '\v':
            case '\t':
            case '\f':
            case '(':
            case ']':
                phy_errno = 2;
                return PHY_ERR;
            case '\0':
                phy_errno = 4;
                return PHY_ERR;
            default:
                if (ctx->n == ctx->nAlloc)
                {
                    ctx->nAlloc += 100;
                    ctx->z = realloc(ctx->z, ctx->nAlloc);
                    if (!ctx->z)
                    {
                        phy_errno = 1;
                        return PHY_ERR;
                    }
                }
                ctx->z[ctx->n++] = c;
        }
    }
    ctx->z[ctx->n] = 0;
    if (ctx->n)
    {
        ctx->p->lab = calloc(ctx->n+1, 1);
        if (!ctx->p->lab)
        {
            phy_errno = 1;
            ctx->p->lab = 0;
            return PHY_ERR;
        }
        strcpy(ctx->p->lab, ctx->z);
    }
    ctx->n = 0;
    return PHY_OK;
}


// Read off a note in a Newick string. Return 0 if successful, 1 on error
static int read_note(struct newick_reader *ctx)
{
    char c;
    int opened;
    c = ctx->newick[ctx->cursor++];
    if (c == '[')
    {
        opened = 1;
        while ((c = ctx->newick[ctx->cursor++]), opened > 0)
        {
            if (c == '\0')
            {
                phy_errno = 4;
                return PHY_ERR;
            }
            if (c == '[')
                ++opened;
            else if (c == ']')
                --opened;
            if (opened)
            {
                if (ctx->n == ctx->nAlloc)
                {
                    ctx->nAlloc += 100;
                    ctx->z = realloc(ctx->z, ctx->nAlloc);
                    if (!ctx->z)
                    {
                        phy_errno = 1;
                        return PHY_ERR;
                    }
                }
                ctx->z[ctx->n++] = c;
            }
        }
        ctx->z[ctx->n] = 0;
        if (ctx->n)
        {
            ctx->p->note = calloc(ctx->n+1, 1);
            if (!ctx->p->note)
            {
                phy_errno = 1;
                ctx->p->note = 0;
                return PHY_ERR;
            }
            strcpy(ctx->p->note, ctx->z);
        }
        ctx->n = 0;
    }
    else
        ctx->cursor--;
    return PHY_OK;
}


// Read off a branch length in a Newick string. Return 0 if successful, 1 on error
static int read_brlen(struct newick_reader *ctx)
{
    char c;
    int toread = 1;
    c = ctx->newick[ctx->cursor++];
    if (c == ':')
    {
        while (toread)
        {
            c = ctx->newick[ctx->cursor++];
            switch (c)
            {
                case 'e':
                case '-':
                case '+':
                case '0':
                case '1':
                case '2':
                case '3':
                case '4':
                case '5':
                case '6':
                case '7':
                case '8':
                case '9':
                case '.':
                    if (ctx->n == ctx->nAlloc)
                    {
                        ctx->nAlloc += 100;
                        ctx->z = realloc(ctx->z, ctx->nAlloc);
                        if (!ctx->z)
                        {
                            phy_errno = 1;
                            return PHY_ERR;
                        }
                    }
                    ctx->z[ctx->n++] = c;
                    break;
                case ',':
                case ')':
                case ';':
                    ctx->cursor--;
                    toread = 0;
                    break;
                case '\0':
                    phy_errno = 4;
                    return PHY_ERR;
                default:
                    phy_errno = 2;
                    return PHY_ERR;
            }
        }
    }
    else
        ctx->cursor--;
    ctx->z[ctx->n] = 0;
    if (ctx->n)
        ctx->p->brlen = atof(ctx->z);
    ctx->n = 0;
    return PHY_OK;
}


// Read a Newick string and return the root node
static struct phy_node *read_newick(struct newick_reader *ctx)
{
    char c;

    if (ctx->newick[strlen(ctx->newick)-1] != ';')
    {
        phy_errno = 4;
        return NULL;
    }

    ctx->nnode++;
    ctx->root = node_new();
    if (ctx->root == NULL)
    {
        phy_errno = 1;
        return NULL;
    }
    ctx->p = ctx->root;
    while ((c = ctx->newick[ctx->cursor++]) != ';')
    {
        switch (c)
        {
            case '(':
                if (ctx->cursor > 1 && !(
                    ctx->newick[ctx->cursor-2] == ','
                    || ctx->newick[ctx->cursor-2] == '('))
                {
                    phy_errno = 4;
                    return NULL;
                }
                ctx->nnode++;
                ctx->q = node_new();
                if (ctx->q == NULL)
                {
                    phy_errno = 1;
                    return NULL;
                }
                phy_node_add_child(ctx->p, ctx->q);
                ctx->p = ctx->q;
                break;
            case ',':
                ctx->nnode++;
                if (!ctx->p->anc && ctx->p != ctx->root)
                {
                    phy_errno = 4;
                    free(ctx->p->lab);
                    free(ctx->p->note);
                    free(ctx->p);
                    return NULL;
                }
                else if (ctx->p == ctx->root)
                {
                    phy_errno = 4;
                    return NULL;
                }
                ctx->q = node_new();
                if (ctx->q == NULL)
                {
                    phy_errno = 1;
                    return NULL;
                }
                phy_node_add_child(ctx->p->anc, ctx->q);
                if (!ctx->p->ndesc)
                    ctx->ntip++;
                ctx->p = ctx->q;
                break;
            case ')':
                if (!ctx->p->ndesc)
                    ctx->ntip++;
                if (!ctx->p->anc && ctx->p != ctx->root)
                {
                    phy_errno = 4;
                    free(ctx->p->lab);
                    free(ctx->p->note);
                    free(ctx->p);
                    return NULL;
                }
                else if (ctx->p == ctx->root)
                {
                    phy_errno = 4;
                    return NULL;
                }
                ctx->p = ctx->p->anc;
                if (ctx->p->ndesc < 2)
                {
                    phy_errno = 3;
                    return NULL;
                }
                break;
            default:
                ctx->cursor--;
                if (read_label(ctx))
                    return NULL;
                if (read_note(ctx))
                    return NULL;
                if (read_brlen(ctx))
                    return NULL;
        }
    }
    if (ctx->p->ndesc < 2)
    {
        phy_errno = 3;
        return NULL;
    }
    return ctx->p;
}


static int write_chars(char *z, struct newick_writer *ctx)
{
    unsigned int len = strlen(z);
    if (len > 0)
    {
        int togo = ctx->nAlloc - ctx->n;
        while (togo <= len)
        {
            ctx->nAlloc += 100;
            ctx->newick = realloc(ctx->newick, ctx->nAlloc);
            if (!ctx->newick)
            {
                phy_errno = 1;
                return PHY_ERR;
            }
            togo = ctx->nAlloc - ctx->n;
        }
        ctx->newick[ctx->n] = 0;
        strcat(ctx->newick, z);
        ctx->n += len;
    }
    return PHY_OK;
}


static int write_brlen(struct phy_node *p, struct newick_writer *ctx)
{
    if (p->brlen > 0)
    {
        char brlen[25] = "";
        sprintf(brlen, ":%f", p->brlen);
        if (write_chars(brlen, ctx))
            return PHY_ERR;
    }
    return PHY_OK;
}


static int write_label(struct phy_node *p, struct newick_writer *ctx)
{
    if (p->lab != 0)
    {
        if (strlen(p->lab) > 0)
        {
            if (write_chars(p->lab, ctx))
                return PHY_ERR;
        }
    }
    return PHY_OK;
}


static int write_note(struct phy_node *p, struct newick_writer *ctx)
{
    if (p->note != 0)
    {
        if (strlen(p->note) > 0)
        {
            if (write_chars(p->note, ctx))
                return PHY_ERR;
        }
    }
    return PHY_OK;
}


static int write_newick(struct phy_node *node, struct newick_writer *ctx)
{
    struct phy_node *d = 0;
    if (node->ndesc)
    {
        for (d = node->lfdesc; d != 0; d = d->next)
        {
            if (d == node->lfdesc)
            {
                if (write_chars("(", ctx))
                    return PHY_ERR;
            }
            write_newick(d, ctx);
            if (d->next)
            {
                if (write_chars(",", ctx))
                    return PHY_ERR;
            }
        }
        if (write_chars(")", ctx))
            return PHY_ERR;
    }
    if (write_label(node, ctx))
        return PHY_ERR;
    if (write_note(node, ctx))
        return PHY_ERR;
    if (write_brlen(node, ctx))
        return PHY_ERR;
    return PHY_OK;
}



/**********************************************************************
**
** Functions accessible via the API
**
**********************************************************************/


int phy_node_alloc(struct phy_node **node)
{
    *node = node_new();
    if (*node)
        return PHY_OK;
    return PHY_ERR;
}


void phy_node_free(struct phy_node *node)
{
    node_free(node);
}


void phy_node_add_child(struct phy_node *parent, struct phy_node *child)
{
    struct phy_node *r;
    switch (parent->ndesc)
    {
        case 0:
            parent->lfdesc = child;
            child->prev = 0;
            break;
        case 1:
            parent->lfdesc->next = child;
            child->prev = parent->lfdesc;
            break;
        case 2:
            parent->lfdesc->next->next = child;
            child->prev = parent->lfdesc->next;
            break;
        case 3:
            parent->lfdesc->next->next->next = child;
            child->prev = parent->lfdesc->next->next;
            break;
        case 4:
            parent->lfdesc->next->next->next->next = child;
            child->prev = parent->lfdesc->next->next->next;
            break;
        case 5:
            parent->lfdesc->next->next->next->next->next = child;
            child->prev = parent->lfdesc->next->next->next->next;
            break;
        default:
            for (r = parent->lfdesc; r->next != 0; r = r->next) {};
            r->next = child;
            child->prev = r;
    }
    parent->ndesc += 1;
    child->anc = parent;
}

// Remove q from p's list of immediate descendants and return q.
// If q is not a child of p return NULL.
struct phy_node *phy_node_prune(struct phy_node *p, struct phy_node *q)
{
    if (!q || !p)
        return NULL;

    if (q->anc != p)
        return NULL;

    struct phy_node *prev = q->prev;
    struct phy_node *next = q->next;

    q->prev = 0;
    q->next = 0;
    q->anc = 0;

    if (prev)
        prev->next = next;
    if (next)
        next->prev = prev;

    if (q == p->lfdesc)
        p->lfdesc = next;

    p->ndesc -=1;

    return q;
}


// Called after a Newick string has been fully processed
struct phy *phy_build(struct phy_node *root, int nnode, int ntip)
{
    int i = 0, j = 0, k = 0;
    struct phy_node *p, *q;
    struct phy *phy = malloc(sizeof(struct phy));
    if (!phy)
    {
        phy_errno = 1;
        return NULL;
    }
    phy->ntip = ntip;
    phy->nnode = nnode;
    phy->root = root;
    phy->nodes = malloc(nnode * sizeof(struct phy_node *));
    if (!phy->nodes)
    {
        phy_errno = 1;
        free(phy);
        return NULL;
    }
    phy->vseq = malloc(nnode * sizeof(int));
    if (!phy->vseq)
    {
        phy_errno = 1;
        free(phy->nodes);
        free(phy);
        return NULL;
    }
    phy->inodes = malloc((nnode-ntip) * sizeof(struct phy_node *));
    if (!phy->inodes)
    {
        phy_errno = 1;
        free(phy->vseq);
        free(phy->nodes);
        free(phy);
        return NULL;
    }

    /* Assign indices to each node and record their visitation sequence
    ** in a preorder traversal beginning from the root. Terminal nodes
    ** are numbered 0 ... ntip-1 and internal nodes are numbered from
    ** ntip ... nnode-1. The root node always has index ntip using this
    ** scheme. */

    p = phy->root;
    while (p)
    {
        phy->nodes[i] = p;
        if (p->ndesc)
        {
            p->index = ntip + j;
            phy->inodes[j++] = p;
        }
        else
        {
            p->index = k++;
        }
        phy->vseq[p->index] = i++;
        if (p->lfdesc)
            p = p->lfdesc;
        else if (p->next)
            p = p->next;
        else
        {
            // on entry p is a terminal node that marks a clade
            // boundary. this will be the last node visited in a
            // preorder traversal of the subtree rooted at each node
            // on the path back from this node up to and including the
            // first encountered node with a ->next sibling
            q = p;
            while (p->anc != 0 && p->next == 0)
            {
                p = p->anc;
                p->lastvisit = q;
            }
            p = p->next;
        }
    }
    return phy;
}


void phy_cursor_prepare_v2(
    struct phy *phy,
    struct phy_node *node,
    struct phy_cursor *cursor,
    int visit,
    int order
){
    cursor->visit = visit;
    cursor->order = order;
    cursor->phy = phy;

    if (order == PREORDER)
    {
        cursor->begin = node;
        cursor->next = node;
        if (node->ndesc)
        {
            if (visit == ALL_NODES)
            {
                cursor->end = phy->nodes[phy->vseq[node->lastvisit->index]];
            }
            else
            {
                struct phy_node *q;
                struct phy_node *p;
                for (p = q = node->lastvisit->anc->lfdesc; q; q = q->next) {
                    if (q->ndesc)
                        p = q;
                }
                while (p->ndesc) {
                    for (p = q = p->lastvisit->anc->lfdesc; q; q = q->next) {
                        if (q->ndesc)
                            p = q;
                    }
                }
                cursor->end = phy->inodes[p->anc->index - phy->ntip];
            }
        }
        else
        {
            cursor->end = node;
        }

        if (node->ndesc)
        {
            if (visit == ALL_NODES)
                cursor->cursor = phy->vseq[node->index] + 1;
            else
                cursor->cursor = node->index - phy->ntip + 1;
        }
    }
    else
    { // POSTORDER
        cursor->end = node;
        if (node->ndesc)
        {
            if (visit == ALL_NODES)
            {
                cursor->begin = phy->nodes[phy->vseq[node->lastvisit->index]];
            }
            else
            {
                struct phy_node *q;
                struct phy_node *p;
                for (p = q = node->lastvisit->anc->lfdesc; q; q = q->next) {
                    if (q->ndesc)
                        p = q;
                }
                while (p->ndesc) {
                    for (p = q = p->lastvisit->anc->lfdesc; q; q = q->next) {
                        if (q->ndesc)
                            p = q;
                    }
                }
                cursor->begin = phy->inodes[p->anc->index - phy->ntip];
            }
        }
        else
        {
            cursor->begin = node;
        }
        cursor->next = cursor->begin;
        if (node->ndesc)
        {
            if (visit == ALL_NODES)
                cursor->cursor = phy->vseq[cursor->begin->index] - 1;
            else
                cursor->cursor = cursor->begin->index - phy->ntip - 1;
        }
    }
}


struct phy_cursor *phy_cursor_prepare(
    struct phy *phy,
    struct phy_node *node,
    int visit,
    int order
){
    phy_cursor_prepare_v2(phy, node, &phy->cursor, visit, order);
    return &phy->cursor;
}


struct phy_node *phy_cursor_step(struct phy_cursor *cursor)
{
    if (!cursor)
        return NULL;

    struct phy_node *node = cursor->next;

    if (node)
    {
        if (node != cursor->end)
        {
            if (cursor->visit == ALL_NODES)
            {
                if (cursor->order == PREORDER)
                    cursor->next = cursor->phy->nodes[cursor->cursor++];
                else
                    cursor->next = cursor->phy->nodes[cursor->cursor--];
            }
            else
            {
                if (cursor->order == PREORDER)
                    cursor->next = cursor->phy->inodes[cursor->cursor++];
                else
                    cursor->next = cursor->phy->inodes[cursor->cursor--];
            }
        }
        else
        {
            cursor->begin = 0;
            cursor->end = 0;
            cursor->next = 0;
        }
    }

    return node;
}


int phy_isbinary(struct phy *phy)
{
    return (phy->nnode ==  2*phy->ntip - 1) ? 1 : 0;
}


int phy_cursor_alloc(struct phy_cursor **cursor)
{
    *cursor = (struct phy_cursor *) malloc(sizeof(struct phy_cursor));

    if (*cursor)
        return PHY_OK;

    phy_errno = 1;
    return PHY_ERR;
}


void phy_cursor_free(struct phy_cursor *cursor)
{
    free(cursor);
}


void phy_free(struct phy *phy)
{
    if (phy)
    {
        int i;
        for (i = 0; i < phy->nnode; ++i)
            node_free(phy->nodes[i]);
        free(phy->nodes);
        free(phy->inodes);
        free(phy->vseq);
        free(phy);
    }
}


struct phy *phy_read_newickstr(const char *newick)
{
    struct phy_node *root = 0;
    struct phy *phy = 0;
    struct newick_reader ctx = {0, 0, 0, 0, 0, 0, newick, 0, 0, 0};
    root = read_newick(&ctx);
    if (root)
    {
        phy = phy_build(root, ctx.nnode, ctx.ntip);
        if (!phy)
            cleanup(ctx.root);
    }
    else
        cleanup(ctx.root);
    free(ctx.z);
    return phy;
}


struct phy *phy_read_newickfile(const char *filename)
{
    struct phy *phy = 0;
    char *newick = 0;
    int nBytes = 0;
    FILE *in = fopen(filename, "r");
    if (in) {
        fseek(in, 0, SEEK_END);
        nBytes = ftell(in);
        rewind(in);
        if (nBytes)
            newick = malloc(nBytes);
        if (newick)
        {
            fread(newick, nBytes, 1, in);
            phy = phy_read_newickstr(newick);
            free(newick);
        }
        fclose(in);
    }
    return phy;
}


char *phy_write_newickstr(struct phy *phy)
{
    struct newick_writer ctx = {0, 0, 0};
    write_newick(phy->root, &ctx);
    write_chars(";", &ctx);
    return ctx.newick;
}


int phy_write_newickfile(
    struct phy *phy, const char *filename, const char *mode)
{
    char *newick = phy_write_newickstr(phy);
    if (!newick)
        return PHY_ERR;
    FILE *out = fopen(filename, mode);
    if (out) {
        fputs(newick, out);
        free(newick);
        fclose(out);
        return PHY_OK;
    }
    return PHY_ERR;
}


struct phy *phy_extract_clade(struct phy_node *node)
{
    struct newick_writer ctx = {0, 0, 0};
    struct phy *phy = 0;
    if (write_newick(node, &ctx))
        return NULL;
    if (write_chars(";", &ctx))
        return NULL;
    phy = phy_read_newickstr(ctx.newick);
    free(ctx.newick);
    if (!phy)
        return NULL;
    phy->root->brlen = 0;
    return phy;
}


struct phy *phy_extract_subtree(
    int ntip,
    struct phy_node **tips,
    struct phy *phy
){
    int i;
    int nnode = 0;
    int bitmask[phy->nnode];
    struct phy_node *p;
    struct phy_node *q;
    struct phy_node *root;
    struct phy_node *head;

    memset(bitmask, 0, phy->nnode * sizeof(int));
    for (i = 0; i < ntip; ++i)
    {
        p = tips[i];
        bitmask[p->index] = 1;
        while ((p = p->anc) != 0)
        {
            if (bitmask[p->index])
                break;
            bitmask[p->index] = 1;
        }
    }

    root = 0;
    head = 0;

    p = phy->root;
    while (p)
    {
        if (bitmask[p->index] && !root)
        {
            nnode++;
            root = node_new();
            if (!root)
            {
                phy_errno = 1;
                return NULL;
            }
            head = root;
        }
        else if (bitmask[p->index])
        {
            nnode++;
            q = node_new();
            if (!q)
            {
                phy_errno = 1;
                cleanup(root);
                return NULL;
            }
            q->brlen = p->brlen;
            if (p->lab)
            {
                q->lab = malloc(strlen(p->lab)+1);
                if (!q->lab)
                {
                    phy_errno = 1;
                    cleanup(root);
                    return NULL;
                }
                strcpy(q->lab, p->lab);
            }
            phy_node_add_child(head, q);
            head = q;
        }
        if (p->lfdesc)
            p = p->lfdesc;
        else if (p->next)
        {
            if (bitmask[p->index])
                head = head->anc;
            p = p->next;
        }
        else
        {
            while (p->anc != 0 && p->next == 0)
            {
                if (bitmask[p->index])
                    head = head->anc;
                p = p->anc;
            }
            if (bitmask[p->index])
                head = head->anc;
            p = p->next;
        }
    }

    root->brlen = 0;
    p = root;
    while (p)
    {
        while (p->ndesc == 1)
        {
            q = p->lfdesc;
            q->brlen += p->brlen;
            q->next = p->next;
            q->prev = p->prev;
            q->anc = p->anc;
            if (p->prev)
                p->prev->next = q;
            if (p->next)
                p->next->prev = q;
            if (p->anc && p->anc->lfdesc == p)
                p->anc->lfdesc = q;
            if (p == root) {
                root = q;
                root->brlen = 0;
            }
            node_free(p);
            p = q;
            nnode--;
        }
        if (p->lfdesc)
            p = p->lfdesc;
        else if (p->next)
            p = p->next;
        else {
            while (p->anc != 0 && p->next == 0)
                p = p->anc;
            p = p->next;
        }
    }

    return phy_build(root, nnode, ntip);
}


void phy_node_rotate(int n, struct phy_node **nodes, struct phy *phy)
{
    int i;
    int j;
    int k;
    struct phy_node *d = 0;
    struct phy_node *p;
    struct phy_node *q;
    struct phy_node *node;

    for (i = 0; i < n; ++i)
    {
        node = nodes[i];

        q = node->lfdesc;
        p = q->next;

        while (q)
        {
            phy_node_prune(node, q);
            if (d)
            {
                q->prev = d;
                d = q;
            }
            else
                d = q;
            q = p;
            if (q)
                p = q->next;
        }

        while (d)
        {
            q = d;
            d = q->prev;
            phy_node_add_child(node, q);
        }

    }

    i = 0;
    j = 0;
    k = 0;

    p = phy->root;
    while (p)
    {
        phy->nodes[i] = p;
        if (p->ndesc)
        {
            p->index = phy->ntip + j;
            phy->inodes[j++] = p;
        }
        else
            p->index = k++;
        phy->vseq[p->index] = i++;

        if (p->lfdesc)
            p = p->lfdesc;
        else if (p->next)
            p = p->next;
        else {
            q = p;
            while (p->anc != 0 && p->next == 0) {
                p = p->anc;
                p->lastvisit = q;
            }
            p = p->next;
        }
    }
}


struct phy *phy_duplicate(struct phy *phy)
{
    char *newick = phy_write_newickstr(phy);
    struct phy *phy_cpy = phy_read_newickstr(newick);
    free(newick);
    return phy_cpy;
}


void phy_reroot(
    struct phy_node *node,
    struct phy **in,
    struct phy **out
){
    int overwrite = *in == *out ? 1 : 0;
    int isrooted = phy_isrooted(*in);
    double b;

    const char *lf_tip;
    const char *rt_tip;

    struct phy *phy;
    struct phy_node *d;
    struct phy_node *p;
    struct phy_node *q = 0;
    struct phy_node *root_on;
    struct phy_node *root;


    phy_node_spanning_pair(node, &lf_tip, &rt_tip);

    phy = phy_duplicate(*in);

    root_on = phy_node_mrca(
        phy, phy_node_find(phy, lf_tip), phy_node_find(phy, rt_tip));

    b = root_on->brlen / 2;
    p = root_on->anc;
    phy_node_prune(p, root_on);
    root = node_new();
    phy_node_add_child(root, root_on);
    root_on->brlen = b;

    // build up a list like so,
    //   * <- 1 <- 2 <- 3 <- current_root
    // where * equals root_on->anc and numbers
    // refer to nodes on the path back to the
    // current root. when re-rooting we need
    // to reverse the direction of each edge on
    // this path.
    while (p->anc)
    {
        d = p->anc;
        // removes p from d's descendants and sets p->prev, p->next
        // and p->anc to NULL
        p = phy_node_prune(d, p);
        if (!q)
            q = p;
        else
        {
            p->prev = q;
            q = p;
        }

        p = d;
    }
    // p is now the current root
    p->prev = q;
    q = p;

    if (isrooted && q->prev)
    {
        // if the tree is currently rooted we want to
        // jump around the current root
        d = q;
        p = q->prev;
        for (d = d->lfdesc; d == p; d = d->next) {}
        phy_node_add_child(p, d);
        d->brlen += p->brlen;
        q = p;
    }
    else if (isrooted)
    {
        // rooting on a node whose parent is the current
        // root, but we still need to jump around the
        // current root. this bit prepares for the
        // statement following the while loop below
        // (which will not be entered if this branch is
        // taken)
        d = q;
        p = root_on;
        for (d = d->lfdesc; d == p; d = d->next) {}
        q = d;
        q->prev = NULL;
        b += d->brlen;
    }

    // reverse the rest of the edges
    while (q->prev) {
        d = q;
        p = q->prev;
        phy_node_add_child(p, d);
        d->brlen = p->brlen;
        q = p;
    }

    phy_node_add_child(root, q);
    q->brlen = b;

    *out = phy_build(root, phy->nnode + (isrooted ? 0 : 1), phy->ntip);

    if (overwrite)
    {
        phy_free(*in);
        *in = *out;
    }

    phy_free(phy);
}


void phy_unroot(struct phy **in, struct phy **out)
{
    int overwrite = *in == *out ? 1 : 0;
    int isrooted = phy_isrooted(*in);

    if (isrooted)
    {
        *out = 0;
        return;
    }

    struct phy *phy;
    struct phy_node *p;
    struct phy_node *q;

    phy = phy_duplicate(*in);

    p = phy_node_prune(phy->root, phy->root->lfdesc);
    q = phy_node_prune(phy->root, phy->root->lfdesc->next);

    q->brlen += p->brlen;
    p->brlen = 0;

    phy_node_add_child(p, q);

    *out = phy_build(p, phy->nnode - 1, phy->ntip);

    if (overwrite)
    {
        phy_free(*in);
        *in = *out;
    }

    phy_free(phy);
}


struct phy_node *phy_node_lfdesc(struct phy_node *node)
{
    return node->lfdesc;
}


struct phy_node *phy_node_rtdesc(struct phy_node *node)
{
    struct phy_node *r;
    switch (node->ndesc)
    {
        case 0:
            return NULL;
        case 2:
            return node->lfdesc->next;
        case 3:
            return node->lfdesc->next->next;
        case 4:
            return node->lfdesc->next->next->next;
        case 5:
            return node->lfdesc->next->next->next->next;
        default:
            for (r = node->lfdesc; r->next != 0; r = r->next) {};
            return r;
    }
}


struct phy_node *phy_node_next(struct phy_node *node)
{
    return node->next;
}


struct phy_node *phy_node_prev(struct phy_node *node)
{
    return node->prev;
}


struct phy_node *phy_node_anc(struct phy_node *node)
{
    return node->anc;
}


double phy_node_brlen(struct phy_node *node)
{
    return node->brlen;
}


int phy_node_index(struct phy_node *node)
{
    return node->index;
}


int phy_node_ndesc(struct phy_node *node)
{
    return node->ndesc;
}


int phy_node_istip(struct phy_node *node)
{
    return node->ndesc ? 0 : 1;
}


int phy_isrooted(struct phy *phy)
{
    return phy->root->ndesc > 2 ? 0 : 1;
}


int phy_nnode(struct phy *phy)
{
    return phy->nnode;
}


int phy_ntip(struct phy *phy)
{
    return phy->ntip;
}


void phy_node_swap(struct phy_node *a, struct phy_node *b)
{
    if (!a->anc || !b->anc || a->anc != b->anc)
        return;

    struct phy_node *prev;
    struct phy_node *next;

    next = a->next;
    prev = a->prev;

    a->next = b->next;
    a->prev = b->prev;
    if (b->next)
        b->next->prev = a;
    if (b->prev)
        b->prev->next = a;
    b->next = next;
    b->prev = prev;
    if (next)
        next->prev = b;
    if (prev)
        prev->next = b;
    if (a->anc->lfdesc == a)
        a->anc->lfdesc = b;
}


void phy_ladderize(struct phy *phy, int *perm)
{
    int i = 0, j = 0, k = 0, sorted;
    struct phy_node *d;
    struct phy_node *p;
    struct phy_node *q;
    struct phy_cursor *cursor;

    cursor = phy_cursor_prepare(
        phy, phy->root, INTERNAL_NODES_ONLY, PREORDER);

    while ((p = phy_cursor_step(cursor)) != 0)
    {

        do {
            sorted = 0;

            d = p->lfdesc;
            while (d && d->next)
            {
                if (d->ndesc > d->next->ndesc)
                {
                    ++sorted;
                    phy_node_swap(d, d->next);
                }
                d = d->next;
            }

        } while (sorted > 0);
    }

    p = phy->root;
    while (p)
    {
        phy->nodes[i] = p;
        if (p->ndesc)
        {
            perm[phy->ntip + j] = p->index;
            p->index = phy->ntip + j;
            phy->inodes[j++] = p;
        }
        else
        {
            perm[k] = p->index;
            p->index = k++;
        }
        phy->vseq[p->index] = i++;

        if (p->lfdesc)
            p = p->lfdesc;
        else if (p->next)
            p = p->next;
        else {
            q = p;
            while (p->anc != 0 && p->next == 0) {
                p = p->anc;
                p->lastvisit = q;
            }
            p = p->next;
        }
    }
}


struct phy_node *phy_root(struct phy *phy)
{
    return phy->root;
}


struct phy_node *phy_node_get(struct phy *phy, int index)
{
    if (index >= 0 && index < phy->nnode)
        return phy->nodes[phy->vseq[index]];
    return NULL;
}


struct phy_node *phy_node_find(struct phy *phy, const char *label)
{
    struct phy_cursor *cursor;
    struct phy_node *node = 0;

    cursor = phy_cursor_prepare(phy, phy->root, ALL_NODES, POSTORDER);
    while ((node = phy_cursor_step(cursor)) != 0)
    {
        if (node->lab)
        {
            if (strcmp(label, node->lab) == 0)
                break;
        }
    }

    return node;
}


void phy_node_add_data(
    struct phy_node *node,
    void *data,
    void (*data_free)(void *)
){
    if (node->data && node->data_free)
        node->data_free(node->data);
    node->data = data;
    node->data_free = data_free;
}


void *phy_node_data(struct phy_node *node)
{
    return node->data;
}


void phy_node_spanning_pair(
    struct phy_node *node,
    const char **a,
    const char **b
){
    if (!node->ndesc)
    {
        *a = node->lab;
        *b = node->lab;
    }
    else
    {
        *b = node->lastvisit->lab;
        node = node->lfdesc;
        while (node->lfdesc)
            node = node->lfdesc;
        *a = node->lab;
    }
}


struct phy_node *phy_node_mrca(
    struct phy *phy,
    struct phy_node *a,
    struct phy_node *b
){
    if (a == b)
        return a;

    int gotcha[phy->nnode];

    memset(gotcha, 0, phy->nnode * sizeof(int));

    while ((a = a->anc) != 0)
        gotcha[a->index] = 1;

    while (b && !gotcha[b->index])
        b = b->anc;

    return b;
}


void phy_node_foreach(
    struct phy *phy,
    struct phy_node *node,
    int visit,
    int order,
    void (*FUN)(struct phy_node *node, struct phy *phy, void *param),
    void *param
){
    struct phy_cursor cursor;
    phy_cursor_prepare_v2(phy, node, &cursor, visit, order);
    while ((node = phy_cursor_step(&cursor)) != 0)
        FUN(node, phy, param);
}


const char *phy_node_label(struct phy_node *node)
{
    return node->lab;
}


const char *phy_node_note(struct phy_node *node)
{
    return node->note;
}


const char *phy_errmsg()
{
    switch (phy_errno)
    {
        case 1:
            phy_errno = 0;
            return PHY_ERR1;
        case 2:
            phy_errno = 0;
            return PHY_ERR2;
        case 3:
            phy_errno = 0;
            return PHY_ERR3;
        case 4:
            phy_errno = 0;
            return PHY_ERR4;
        default:;
    }
    return "no errors detected";
}
