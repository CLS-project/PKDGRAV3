#ifndef RBTREE_H
#define RBTREE_H

typedef struct rb_node {
    struct rb_node *link[2];
    int red;
    } RB_NODE;
typedef RB_NODE *RB_TREE;

typedef struct {
    void *ctx;
    int (*cmp)(void *ctx, const void *a, const void *b);
    RB_NODE *(*newnode)(void *ctx, const void *data);
    void (*delnode)(void *ctx, RB_NODE *node);
    int iSize;
    } RB_TYPE;

void rb_type_create(
    RB_TYPE *rbt,int iSize, void *ctx,
    int (*cmp)(void*,const void *,const void *),/* Required */
    RB_NODE *(*newnode)(void*,const void*),     /* Defaults to malloc */
    void (*delnode)(void*,RB_NODE*) );          /* Defaults to free */

void rb_free( const RB_TYPE *rbt, RB_TREE *root );

int rb_size( RB_TREE *root );

RB_NODE *rb_search(const RB_TYPE *rbt, RB_NODE *n, void *data);
int rb_insert( const RB_TYPE *rbt, RB_TREE *root, void *data );
int rb_remove ( const RB_TYPE *rbt, RB_TREE *root, void *data );
#endif
