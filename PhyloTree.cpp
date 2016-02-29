#include <iostream>
#include <fstream>
#include <stdlib.h>
#include "PhyloTree.h"
#include "Exceptions.h"
#include "VerbosityLevels.h"

PhyloTree::~PhyloTree() {
  mTreeRoot.clearNode();
  mLeavesSpecies.clear();
  mInternalNodes.clear();
}

void PhyloTree::fillSpecies(TreeNode *aNode) {
  if (aNode->isLeaf()) {
    mLeavesSpecies.push_back(aNode);
  } else {
    TreeNode *m;
    for (unsigned int idx = 0; (m = aNode->getChild(idx)) != NULL; ++idx)
      fillSpecies(m);
  }
}

const std::vector<std::string> &PhyloTree::getSpecies(void) const {
  mSpeciesList.clear();

  std::vector<TreeNode *>::const_iterator is(mLeavesSpecies.begin());
  const std::vector<TreeNode *>::const_iterator end(mLeavesSpecies.end());
  for (; is != end; ++is) {
    std::string label = (*is)->getLabel();
    mSpeciesList.push_back(label);
  }

  return mSpeciesList;
}

void PhyloTree::fillInternalBranches(TreeNode *aNode) {
  if (/*!aNode->isLeaf() && */ aNode != &mTreeRoot) {
    mInternalNodes.push_back(aNode);
  }

  TreeNode *m;
  for (unsigned int idx = 0; (m = aNode->getChild(idx)) != NULL; ++idx)
    fillInternalBranches(m);
}

size_t PhyloTree::getMarkedInternalBranch(void) const {
  const size_t nin = mInternalNodes.size();

  // std::cout << "NIN " << nin << std::endl;

  size_t marked_branch = 0;
  for (; marked_branch < nin; ++marked_branch) {
    if (!mInternalNodes[marked_branch]->getType().empty())
      break;
  }

  if (marked_branch >= nin)
    return UINT_MAX;
  return marked_branch;
}

/*std::set<size_t> PhyloTree::getMarkedBranches(void) const
 {
 const size_t nin = mInternalNodes.size();

 std::cout << "NIN " << nin << std::endl;

 std::set<size_t> marked_branches;
 size_t current_branch = 0;
 for(; current_branch < nin; ++current_branch)
 {
 if(!mInternalNodes[current_branch]->getType().empty())
 marked_branches.insert(current_branch);
 }

 return marked_branches;
 }*/

unsigned int PhyloTree::cloneTree(ForestNode *aForestNode, unsigned int aTreeId,
                                  size_t aNumSites,
                                  CacheAlignedDoubleVector &aProbVectors,
                                  const TreeNode *aTreeNode,
                                  unsigned int aNodeId) const {
  unsigned int id;

  // Start with the tree root
  if (aTreeNode == NULL) {
    aTreeNode = &mTreeRoot;
    aNodeId = UINT_MAX;
    aForestNode->mParent = NULL;
    id = 0;
  } else {
    id = aNodeId + 1;
  }

  // Set the root values
  aForestNode->mBranchId = aNodeId;
  aForestNode->mOwnTree = aTreeId;

  // Set the internal branch identifier
  const size_t nn = mInternalNodes.size();
  size_t int_id = 0;
  for (; int_id < nn; ++int_id)
    if (aTreeNode == mInternalNodes[int_id])
      break;
  aForestNode->mInternalNodeId =
      (int_id < nn) ? static_cast<unsigned int>(int_id) : UINT_MAX;

#ifndef NEW_LIKELIHOOD
  // Set the pointers.		The sequence is: Branch -> Set -> Site -> 1:N
  // Set the pointers (best). The sequence is: Branch -> Site -> Set -> 1:N
  for (int i = 0; i < Nt; ++i) {
    // aForestNode->mProb[i] =
    // &aProbVectors[VECTOR_SLOT*(aNumSites*Nt*id+aNumSites*i+aTreeId)];
    aForestNode->mProb[i] =
        &aProbVectors[VECTOR_SLOT * (aNumSites * Nt * id + aTreeId * Nt + i)];
  }
#endif

  // Recurse
  TreeNode *m = NULL;
  for (unsigned int idx = 0; (m = aTreeNode->getChild(idx)) != NULL; ++idx) {
    ForestNode *rn = new ForestNode;
    aForestNode->mChildrenList.push_back(rn);
    aForestNode->mChildrenCount++;
#ifndef NEW_LIKELIHOOD
    aForestNode->mOtherTreeProb.push_back(NULL);
#endif
    rn->mParent = aForestNode;
    id = cloneTree(rn, aTreeId, aNumSites, aProbVectors, m, id);
  }

  return id;
}

unsigned int PhyloTree::collectGlobalTreeData(
    std::vector<std::string> &aNodeNames, std::vector<double> &aBranchLengths,
    std::set<int> &aInternalBranches, size_t *aMarkedIntBranch,
    std::set<int> &aMarkedBranches, const TreeNode *aTreeNode,
    unsigned int aNodeId) const {
  unsigned int id;

  // Start with the tree root
  if (aTreeNode == NULL) {
    aTreeNode = &mTreeRoot;
    aNodeId = UINT_MAX;
    id = 0;
    *aMarkedIntBranch = getMarkedInternalBranch();
    // aMarkedBranches.insert(*aMarkedIntBranch); // Get all marked branches
  } else {
    id = aNodeId + 1;
  }

  // Save the node values
  if (!aTreeNode->getType().empty())
    aMarkedBranches.insert(aNodeId); // Get all marked branches
  aNodeNames.push_back(aTreeNode->getLabel());
  aBranchLengths.push_back(aTreeNode->getLen());
  if (!aTreeNode->isLeaf())
    aInternalBranches.insert(aNodeId);

  // Recurse
  TreeNode *m;
  for (unsigned int idx = 0; (m = aTreeNode->getChild(idx)) != NULL; ++idx) {
    id = collectGlobalTreeData(aNodeNames, aBranchLengths, aInternalBranches,
                               aMarkedIntBranch, aMarkedBranches, m, id);
  }

  return id;
}

void PhyloTree::countNullBranchLengths(int &aOnLeafCnt, int &aOnIntCnt,
                                       const TreeNode *aTreeNode) const {
  // Start with the tree root (don't check its branch length that is invalid)
  if (aTreeNode == NULL) {
    aTreeNode = &mTreeRoot;
  } else {
    // Check branch length. It should not be null for leaves
    if (aTreeNode->getLen() < 1e-14) {
      if (aTreeNode->isLeaf())
        ++aOnLeafCnt;
      else
        ++aOnIntCnt;
    }
  }

  // Recurse
  TreeNode *m;
  for (unsigned int idx = 0; (m = aTreeNode->getChild(idx)) != NULL; ++idx) {
    countNullBranchLengths(aOnLeafCnt, aOnIntCnt, m);
  }
}

void PhyloTree::checkRootBranches(void) // const
{
  TreeNode *m;
  unsigned int cnt_root_branches = 0;
  unsigned int cnt_root_leaves = 0;
  for (; (m = mTreeRoot.getChild(cnt_root_branches)) != NULL;
       ++cnt_root_branches) {
    if (m->isLeaf())
      ++cnt_root_leaves;
  }

  if (mVerboseLevel >= VERBOSE_INFO_OUTPUT) {
    std::cout << std::endl
              << "Root has " << cnt_root_branches << " children of which "
              << cnt_root_leaves << " are leaves" << std::endl;
  }

  if (cnt_root_branches == 2) {
    if (mVerboseLevel >= VERBOSE_ONLY_RESULTS)
      std::cout << std::endl
                << "This is a rooted tree. Tree is unrooted !" << std::endl;

    TreeNode *c0, *c1, *cx; // cx is the child with at least 2 children
    c0 = mTreeRoot.getChild(0);
    c1 = mTreeRoot.getChild(1);

    // get the length and marking of new branch

    double new_len = c0->getLen() + c1->getLen();

    // std::cout << " Length Root : " << mTreeRoot.getLen() << std::endl;
    // std::cout << " Length Child 0 : " << c0->getLen() << std::endl;
    // std::cout << " Length Child 1 : " << c1->getLen() << std::endl;
    // std::cout << " New Length : " << new_len << std::endl;

    // bool marked=0;
    // if (c0->getType()!=NULL || c1->getType()!=NULL ) marked = true;

    // choose the child who has also two children and adjust lengths and marks

    // std::cout << " type of c0 : " << c0->getType() << std::endl;
    // std::cout << " type of c1 : " << c1->getType() << std::endl;

    std::string newBranchType;

    if (c0->getType().size() > 0)
      newBranchType = c0->getType();
    else
      newBranchType = c1->getType();

    if (c0->getChild(0) != NULL && c0->getChild(1) != NULL)

    {
      cx = c0;
      c1->addLen(new_len);
      c1->addType(newBranchType);
      mTreeRoot.delChild(0);
    }

    else if (c1->getChild(0) != NULL && c1->getChild(1) != NULL)

    {
      cx = c1;
      c0->addLen(new_len);
      c0->addType(newBranchType);
      mTreeRoot.delChild(1);
    } else

      throw FastCodeMLFatal(
          "Both children of root have only one child. Invalid tree. Quitting.");

    // std::cout << " len of cx : " << cx->getLen() << std::endl;

    // remove cx from stack of children of root
    // all children of cx will be children of root
    // parent of all children of cx will become root

    TreeNode *t;
    unsigned int cx_children = 0;

    for (; (t = cx->getChild(cx_children)) != NULL; cx_children++) {
      mTreeRoot.addChild(t);
      t->addParent(&mTreeRoot);
    }

    // refill the tree

    // Fill the leaves
    mLeavesSpecies.clear();
    this->fillSpecies(&mTreeRoot);

    // Fill internal branches
    mInternalNodes.clear();
    this->fillInternalBranches(&mTreeRoot);

    // cx->clearNode();

    // c0->addLen(new_len);

    // c1->addParent(c0);
    // c1->addLen(new_len);

    // mTreeRoot=*c0;

    // this->mTreeRoot=*c0;
    // mTreeRoot.clearNode();

    // mTreeRoot.clearNode();

    // child one becomes child of child zero

    // adjust the new length and mark of the new branch between child zero and
    // child one

    // TreeNode *m;
    // unsigned int cnt_root_branches = 0;
    // unsigned int cnt_root_leaves   = 0;
    // for(; (m = mTreeRoot.getChild(cnt_root_branches)) != NULL;
    // ++cnt_root_branches)
    //{
    //	if(m->isLeaf()) ++cnt_root_leaves;
    //}

    // std::cout << std::endl << "Root has " << cnt_root_branches << " children
    // of which " << cnt_root_leaves << " are leaves" << std::endl;
  }
  // if it is an invalid tree then raise exception
  if (cnt_root_branches < 2)
    throw FastCodeMLFatal("Root has only one branch. Invalid tree. Quitting.");

  // if it is an invalid tree then raise exception
  if (cnt_root_branches == cnt_root_leaves)
    throw FastCodeMLFatal(
        "Root points only to leaves. Invalid tree. Quitting.");
}
