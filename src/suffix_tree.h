/******************************************************************************
Suffix Tree Version 2.1
by:         Dotan Tsadok.
Instructor: Mr. Shlomo Yona.
University of Haifa, Israel.
December 2002.

Current maintainer:
	Shlomo Yona	<shlomo@cs.haifa.ac.il>

DESCRIPTION OF THIS FILE:
This is the decleration file suffix_tree.h and it contains declarations of the
interface functions for constructing, searching and deleting a suffix tree, and
to data structures for describing the tree.

COPYRIGHT
Copyright 2002-2003 Shlomo Yona

LICENSE
This library is free software; you can redistribute it and/or modify it
under the same terms as Perl itself.
*******************************************************************************/

#if !defined(SUFFIXTREE_H)
#define SUFFIXTREE_H

/* A type definition for a 32 bits variable - a double word. */
#define     DBL_WORD      unsigned long   

/* Error return value for some functions. Initialized  in ST_CreateTree. */
DBL_WORD    ST_ERROR;

/******************************************************************************/
/*                           DATA STRUCTURES                                  */
/******************************************************************************/
/* This structure describes a node and its incoming edge */
typedef struct SUFFIXTREENODE
{
   /* A linked list of sons of that node */
   struct SUFFIXTREENODE*   sons;
   /* A linked list of right siblings of that node */
   struct SUFFIXTREENODE*   right_sibling;
   /* A linked list of left siblings of that node */
   struct SUFFIXTREENODE*   left_sibling;
   /* A pointer to that node's father */
   struct SUFFIXTREENODE*   father;
   /* A pointer to the node that represents the largest 
   suffix of the current node */
   struct SUFFIXTREENODE*   suffix_link;
   /* Index of the start position of the node's path */
   DBL_WORD                 path_position;
   /* Start index of the incoming edge */
   DBL_WORD                 edge_label_start;
   /* End index of the incoming edge */
   DBL_WORD                 edge_label_end;
   /* A counter for search hits on this node */
   /* (C) Jens Kleinjung, London 2006) */
	int hit;
	float entropy;
	_Bool symbol; /* if hit, this symbol is part of the alphabet */
} NODE;

/* This structure describes a suffix tree */
typedef struct SUFFIXTREE
{
   /* The virtual end of all leaves */
   DBL_WORD                 e;
   /* The one and only real source string of the tree. All edge-labels
      contain only indices to this string and do not contain the characters
      themselves */
   char*           tree_string;
   /* The length of the source string */
   DBL_WORD                 length;
   /* The node that is the head of all others. It has no siblings nor a
      father */
   NODE*                    root;
   /* A counter for search hit/miss on this tree */
   /* (C) Jens Kleinjung, London 2006) */
   int allhit, allmiss;
	float allentropy; /* summed entropy over all nodes */
    int allsymbol; /* number of non-zero hit nodes =
                    number of symbols in alphabet to normalise entropy */
} SUFFIX_TREE;


/******************************************************************************/
/*                         INTERFACE FUNCTIONS                                */
/******************************************************************************/
/* 
   ST_CreateTree :
   Allocates memory for the tree and starts Ukkonen's construction algorithm by
   calling SPA n times, where n is the length of the source string.

   Input : The source string and its length. The string is a sequence of
           characters (maximum of 256 different symbols) and not
           null-terminated. The only symbol that must not appear in the string
           is $ (the dollar sign). It is used as a unique symbol by the
           algorithm and is appended automatically at the end of the string (by
           the program, not by the user!). The meaning of the $ sign is
           connected to the implicit/explicit suffix tree transformation,
           detailed in Ukkonen's algorithm.

   Output: A pointer to the newly created tree. Keep this pointer in order to
           perform operations like search and delete on that tree. Obviously,
           no de-allocating of the tree space could be done if this pointer is
           lost, as the tree is allocated dynamically on the heap.
*/

SUFFIX_TREE* ST_CreateTree(const char*   str, DBL_WORD length);

/******************************************************************************/
/*
   ST_FindSubstring :
   Traces for a string in the tree. This function is used for substring search
   after tree construction is done. It simply traverses down the tree starting
   from the root until either the searched string is fully found ot one
   non-matching character is found. In this function skipping is not enabled
   because we don't know wether the string is in the tree or not (see function
   trace_string above).

   Input : The tree, the string W, and the length of W.

   Output: If the substring is found - returns the index of the starting
           position of the substring in the tree source string. If the substring
           is not found - returns ST_ERROR.
*/

DBL_WORD ST_FindSubstring(SUFFIX_TREE*      tree,   /* The suffix array */
                          char*    W,      /* The substring to find */
                          DBL_WORD          P);     /* The length of W */

/******************************************************************************/
/*
   ST_PrintTree :
   This function prints the tree. It simply starts the recoursive function
   ST_PrintNode from the root

   Input : The tree to be printed.
  
   Output: A print out of the tree to the screen.
*/

void ST_PrintTree(SUFFIX_TREE* tree);

/******************************************************************************/
/*
   ST_DeleteTree
   Deletes a whole suffix tree by starting a recoursive call to ST_DeleteSubTree
   from the root. After all of the nodes have been deleted, the function deletes
   the structure that represents the tree.

   Input : The tree to be deleted.

   Output: None.
*/

void ST_DeleteTree(SUFFIX_TREE* tree);

/******************************************************************************/
/*
   ST_SelfTest
   Self test of the tree - search for all substrings of the main string. See
   testing paragraph in the readme.txt file.

   Input : The tree to test.

   Output: 1 for success and 0 for failure. Prints a result message to the
           screen.
*/

DBL_WORD ST_SelfTest(SUFFIX_TREE* tree);

/******************************************************************************/
/*
	ST_InitNodeHits :
	This function initialies all hit counters in the tree.
	Derived from 'ST_PrintTree' function.
	(C) Jens Kleinjung, London 2006

   Input : The tree, the node that is the root of the subtree, and the depth of 
           that node. The node itself is initialised.
           In each recursive call, the depth is increased.
  
   Output: No output.
*/

void ST_InitNodeHits(SUFFIX_TREE* tree, NODE* node1, long depth);

/******************************************************************************/
/*
	ST_InitHitsTree :
	This function initialies all hit counters in the tree.
	Derived from 'ST_PrintTree' function.
	(C) Jens Kleinjung, London 2006

   Input : The tree to be initialised.
  
   Output: No output.
*/

void ST_InitTreeHits(SUFFIX_TREE* tree);

/******************************************************************************/
/*
	ST_NodeEntropy :
	This function calculates substring entropy based on the node hits.
	Derived from 'ST_PrintTree' function.
	(C) Jens Kleinjung, London 2006

   Input : The tree, the node that is the root of the subtree, and the depth of 
           that node. The hit number of the node itself is used.
           In each recursive call, the depth is increased.
  
   Output: No output.
*/

void ST_NodeEntropy(SUFFIX_TREE* tree, NODE* node1, long depth);

/******************************************************************************/
/*
	ST_TreeEntropy :
	This function calculates the substring entropy of the tree based on all node hits.
	Derived from 'ST_PrintTree' function.
	(C) Jens Kleinjung, London 2006

   Input : The tree to process.
  
   Output: No output.
*/

float ST_TreeEntropy(SUFFIX_TREE* tree);

#endif

