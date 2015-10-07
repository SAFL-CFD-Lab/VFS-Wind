/*****************************************************************
* Copyright (C) by Regents of the University of Minnesota.       *
*                                                                *
* This Software is released under GNU General Public License 2.0 *
* http://www.gnu.org/licenses/gpl-2.0.html                       *
*                                                                *
******************************************************************/

typedef struct node {
  PetscInt Node;
  struct node *next;
} node;

typedef struct list{
  node *head;
} List;

void initlist(List *ilist) {
  ilist->head = PETSC_NULL;
}

void insertnode(List *ilist, PetscInt Node)
{
  node *new;
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
    PetscMalloc(sizeof(node), &new);
    new->next = ilist->head;
    new->Node = Node;
    ilist->head = new;
  }
}

void destroy(List *ilist)
{
  node *current;
  while (ilist->head) {
    current = ilist->head->next;
    PetscFree(ilist->head);
    ilist->head = current;
  }
}


typedef struct list_node {
  PetscInt	index;
  struct list_node *next;
} Node_List;

typedef struct IBMListNode {
  IBMInfo ibm_intp;
  struct IBMListNode* next;
} IBMListNode;

typedef struct IBMList {
  IBMListNode *head;
} IBMList;

void InitIBMList(IBMList *ilist) {
  ilist->head = PETSC_NULL;
}

void AddIBMNode(IBMList *ilist, IBMInfo ibm_intp)
{
  IBMListNode *new;
  PetscNew(IBMListNode, &new);
  new->next = ilist->head;
  new->ibm_intp = ibm_intp;
  ilist->head = new;
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
