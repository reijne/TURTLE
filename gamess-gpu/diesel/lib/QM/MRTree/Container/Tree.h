//***********************************************************************
//
//	Name:			Tree.h
//
//	Description:	a tree base class
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			14.03.1997
//
//
//
//
//
//***********************************************************************

#ifndef __TREE_H
#define __TREE_H

#include "../../../../config.h"

#include <iostream>
using std::cout;

template <class ParentType>
class Tree {
public:
	Tree() { cout << "Error: default constructor Tree::Tree() called"; exit(1); }
	Tree(const ParentType *parent);
	
	const ParentType * getParent() const;
	void	setParent(const ParentType *);

protected:
const ParentType *	parent;

private:
};

template <class ParentType>
inline
Tree<ParentType>::Tree(const ParentType *_parent)
{	parent = _parent;	}


template <class ParentType>
inline
const ParentType *	Tree<ParentType>::getParent() const
{	return parent;	}

template <class ParentType>
inline
void	Tree<ParentType>::setParent(const ParentType *p)
{	parent = p;	}


#endif
