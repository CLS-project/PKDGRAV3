/*
**  Implementation of Red Black Trees.
**
**  Original code by Julienne Walker who graciously placed it
**  in the public domain.  Tutoral and code can be found at:
**
**  http://eternallyconfuzzled.com/tuts/datastructures/jsw_tut_rbtree.aspx
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "rbtree.h"

#define is_red(node) (node != NULL && node->red != 0)

static RB_NODE *rb_single ( RB_NODE *node, int dir ) {
    RB_NODE *save = node->link[!dir];
    
    node->link[!dir] = save->link[dir];
    save->link[dir] = node;

    node->red = 1;
    save->red = 0;
 
    return save;
    }
 
static RB_NODE *rb_double ( RB_NODE *node, int dir ) {
    node->link[!dir] = rb_single ( node->link[!dir], !dir );
    return rb_single ( node, dir );
    }

void rb_type_create(
    RB_TYPE *rbt, int iSize, void *ctx,
    int (*cmp)(void*,const void *,const void *),
    RB_NODE *(*newnode)(void*,const void*),
    void (*delnode)(void*,RB_NODE*) ) {
    rbt->ctx = ctx;
    rbt->cmp = cmp;
    rbt->newnode = newnode;
    rbt->delnode = delnode;
    rbt->iSize = iSize;
    assert(cmp != NULL);
    assert(iSize > sizeof(RB_NODE));
    }

void rb_free( const RB_TYPE *rbt, RB_TREE *root ) {
    RB_NODE *node = *root;
    if ( node != NULL ) {
	rb_free(rbt, &node->link[0]);
	rb_free(rbt, &node->link[1]);
	if ( !rbt->delnode ) free(node);
	else (*rbt->delnode)(rbt->ctx,node);
	*root = NULL;
	}
    }

void docount( RB_NODE *node, int *piSize ) {
    if ( node != NULL ) {
	(*piSize)++;
	docount(node->link[0],piSize);
	docount(node->link[1],piSize);
	}
    }


int rb_size( RB_TREE *root ) {
    int iSize = 0;
    docount(*root,&iSize);
    return iSize;
    }

RB_NODE *rb_search(const RB_TYPE *rbt, RB_NODE *n, void *data) {
    int cmp, dir;
    while(n!=NULL) {
	cmp = (*rbt->cmp)(rbt->ctx,n+1,data);
	if ( cmp == 0 ) break;
	dir = cmp < 0;
	n = n->link[dir];
	}
    return n;
    }

static RB_NODE *make_node(const RB_TYPE *rbt, void *data) {
    RB_NODE *n;
    if ( !rbt->newnode) {
	n = malloc(rbt->iSize);
	memcpy(n+1,data,rbt->iSize-sizeof(RB_NODE));
	}
    else n = (*rbt->newnode)(rbt->ctx,data);
    //tree->nNodes++;
    n->red = 1; /* 1 is red, 0 is black */
    n->link[0] = NULL;
    n->link[1] = NULL;
    return n;
    }

int rb_insert( const RB_TYPE *rbt, RB_TREE *root, void *data ) {
    RB_NODE *node = *root;
    int cmp;
    int ins = 0;
    if ( node == NULL ) {
	/* Empty tree case */
	*root = node = make_node(rbt,data);
	assert(node!=NULL);
	ins = 1;
	}
    else {
	RB_NODE head = {{0,0},0}; /* False tree root */
 
	RB_NODE *g, *t;     /* Grandparent & parent */
	RB_NODE *p, *q;     /* Iterator & parent */
	int dir = 0, last=0;
 
	/* Set up helpers */
	t = &head;
	g = p = NULL;
	q = t->link[1] = node;
 
	/* Search down the tree */
	for ( ; ; ) {
	    if ( q == NULL ) {
		/* Insert new node at the bottom */
		p->link[dir] = q = make_node(rbt,data);
		assert(q!=NULL);
		ins = 1;
		}
	    else if ( is_red ( q->link[0] ) && is_red ( q->link[1] ) ) {
		/* Color flip */
		q->red = 1;
		q->link[0]->red = 0;
		q->link[1]->red = 0;
		}
 
	    /* Fix red violation */
	    if ( is_red ( q ) && is_red ( p ) ) {
		int dir2 = t->link[1] == g;
 
		if ( q == p->link[last] )
		    t->link[dir2] = rb_single ( g, !last );
		else
		    t->link[dir2] = rb_double ( g, !last );
		}
 
	    /* Stop if found */
	    cmp = (*rbt->cmp)(rbt->ctx,q+1,data);
	    if ( cmp == 0 )
		break;
 
	    last = dir;
	    dir = cmp < 0;
 
	    /* Update helpers */
	    if ( g != NULL )
		t = g;
	    g = p, p = q;
	    q = q->link[dir];
	    }
 
	/* Update root */
	*root = node = head.link[1];
	}
 
    /* Make root black */
    node->red = 0;
 
    return ins;
    }

int rb_remove ( const RB_TYPE *rbt, RB_TREE *root, void *data ) {
    RB_NODE *node = *root;
    if ( node != NULL ) {
	RB_NODE head = {{0,0},0}; /* False tree root */
	RB_NODE *q, *p, *g; /* Helpers */
	RB_NODE *f = NULL;  /* Found item */
	int dir = 1, cmp;
 
	/* Set up helpers */
	q = &head;
	g = p = NULL;
	q->link[1] = node;
 
	/* Search and push a red down */
	while ( q->link[dir] != NULL ) {
	    int last = dir;
 
	    /* Update helpers */
	    g = p, p = q;
	    q = q->link[dir];
	    cmp = (*rbt->cmp)(rbt->ctx,q+1,data);
	    dir = cmp < 0;
 
	    /* Save found node */
	    if ( cmp == 0 )
		f = q;
 
	    /* Push the red node down */
	    if ( !is_red ( q ) && !is_red ( q->link[dir] ) ) {
		if ( is_red ( q->link[!dir] ) )
		    p = p->link[last] = rb_single ( q, dir );
		else if ( !is_red ( q->link[!dir] ) ) {
		    RB_NODE *s = p->link[!last];
 
		    if ( s != NULL ) {
			if ( !is_red ( s->link[!last] ) && !is_red ( s->link[last] ) ) {
			    /* Color flip */
			    p->red = 0;
			    s->red = 1;
			    q->red = 1;
			    }
			else {
			    int dir2 = g->link[1] == p;
 
			    if ( is_red ( s->link[last] ) )
				g->link[dir2] = rb_double ( p, last );
			    else if ( is_red ( s->link[!last] ) )
				g->link[dir2] = rb_single ( p, last );
 
			    /* Ensure correct coloring */
			    q->red = g->link[dir2]->red = 1;
			    g->link[dir2]->link[0]->red = 0;
			    g->link[dir2]->link[1]->red = 0;
			    }
			}
		    }
		}
	    }
 
	/* Replace and remove if found */
	if ( f != NULL ) {
	    if (f!=q) memcpy(f+1,q+1,rbt->iSize-sizeof(RB_NODE));
	    p->link[p->link[1] == q] = q->link[q->link[0] == NULL];
	    if ( !rbt->delnode ) free(q);
	    else (*rbt->delnode)(rbt->ctx,q);
	    }
 
	/* Update root and make it black */
	*root = node = head.link[1];
	if ( node != NULL )
	    node->red = 0;
	}
 
    return 1;
    }

int rb_assert ( RB_TYPE *rbt, RB_NODE *node ) {
    int lh, rh;
 
    if ( node == NULL )
	return 1;
    else {
	RB_NODE *ln = node->link[0];
	RB_NODE *rn = node->link[1];
 
	/* Consecutive red links */
	if ( is_red ( node ) ) {
	    if ( is_red ( ln ) || is_red ( rn ) ) {
		puts ( "Red violation" );
		return 0;
		}
	    }
 
	lh = rb_assert ( rbt, ln );
	rh = rb_assert ( rbt, rn );
 
	/* Invalid binary search tree */
	if ( ( ln != NULL && (*rbt->cmp)(rbt->ctx,ln+1,node+1) >= 0 )
	     || ( rn != NULL && (*rbt->cmp)(rbt->ctx,rn+1,node+1) <= 0 ))
	    {
	    puts ( "Binary tree violation" );
	    return 0;
	    }
 
	/* Black height mismatch */
	if ( lh != 0 && rh != 0 && lh != rh ) {
	    puts ( "Black violation" );
	    return 0;
	    }
 
	/* Only count black links */
	if ( lh != 0 && rh != 0 )
	    return is_red ( node ) ? lh : lh + 1;
	else
	    return 0;
	}
    }

#ifdef TEST_RBTREE
typedef struct {
    RB_NODE node;
    int value;
    } VNODE;

int cmp(void *ctx, const void *avoid,const void *bvoid) {
    const int *a = (int *)(avoid);
    const int *b = (int *)(bvoid);
    return *a - *b;
    }

void print_tree(RB_NODE *node) {
    VNODE *v = (VNODE *)(node);
    if ( node == NULL ) return;
    print_tree(node->link[0]);
    printf( "%d\n", v->value);
    print_tree(node->link[1]);
    }

int main() {
    RB_TREE root = NULL;
    RB_TYPE rbt;
    int v,i;

    rb_type_create( &rbt, sizeof(VNODE), 0, cmp, NULL, NULL );

    srand(2);
    for( i=0; i<1000; i++ ) {
	v = rand();
	//v = 1000000-i;
	rb_insert( &rbt, &root, &v );
	}
    rb_assert(&rbt,root);

    printf( "Tree has %d elements\n", rb_size(&root) );

    srand(2);
    for( i=0; i<1000; i++ ) {
	v = rand();
	//v = 1000000-i;
	rb_remove( &rbt, &root, &v );
	}
    rb_assert(&rbt,root);
    printf( "Tree has %d elements\n", rb_size(&root) );

    rb_free(&rbt,&root);

    return 0;
    }
#endif
