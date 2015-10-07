/*****************************************************************
* Copyright (C) by Regents of the University of Minnesota.       *
*                                                                *
* This Software is released under GNU General Public License 2.0 *
* http://www.gnu.org/licenses/gpl-2.0.html                       *
*                                                                *
******************************************************************/

#include "variables.h"

void initlist(LIST *ilist) {
  ilist->head = PETSC_NULL;
}

void insertnode(LIST *ilist, PetscInt Node)
{
  node *_new;
  node *current;
  current = ilist->head;

  PetscTruth Exist = PETSC_FALSE;
  while(current) {
    if (Node == current->Node) {
      Exist = PETSC_TRUE;
    }
    if (Exist) break;
    current = current->next;
  }
  if (!Exist) {
    PetscMalloc(sizeof(node), &_new);
    _new->next = ilist->head;
    _new->Node = Node;
    ilist->head = _new;
  }
}

void destroy(LIST *ilist)
{
  node *current;
  while (ilist->head) {
    current = ilist->head->next;
    PetscFree(ilist->head);
    ilist->head = current;
  }
}

void InitIBMList(IBMList *ilist) {
  ilist->head = PETSC_NULL;
}

void AddIBMNode(IBMList *ilist, IBMInfo ibm_intp)
{
  IBMListNode *_new;
  //Petscnew(IBMListNode, &_new);
	PetscMalloc( sizeof(IBMListNode), &_new);	// seokkoo
  _new->next = ilist->head;
  _new->ibm_intp = ibm_intp;
  ilist->head = _new;
}

void DestroyIBMList(IBMList *ilist)
{
  IBMListNode *current;
  while (ilist->head) {
    current = ilist->head->next;
    PetscFree(ilist->head);
    ilist->head = current;
  }
}

