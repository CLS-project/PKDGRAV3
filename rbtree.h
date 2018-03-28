/*  This file is part of PKDGRAV3 (http://www.pkdgrav.org/).
 *  Copyright (c) 2001-2018 Joachim Stadel & Douglas Potter
 *
 *  PKDGRAV3 is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  PKDGRAV3 is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with PKDGRAV3.  If not, see <http://www.gnu.org/licenses/>.
 */

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
    size_t iSize;
    } RB_TYPE;

void rb_type_create(
    RB_TYPE *rbt,size_t iSize, void *ctx,
    int (*cmp)(void*,const void *,const void *),/* Required */
    RB_NODE *(*newnode)(void*,const void*),     /* Defaults to malloc */
    void (*delnode)(void*,RB_NODE*) );          /* Defaults to free */

void rb_free( const RB_TYPE *rbt, RB_TREE *root );

size_t rb_size( RB_TREE *root );

RB_NODE *rb_search(const RB_TYPE *rbt, RB_NODE *n, void *data);
int rb_insert( const RB_TYPE *rbt, RB_TREE *root, void *data );
int rb_remove ( const RB_TYPE *rbt, RB_TREE *root, void *data );
#endif
