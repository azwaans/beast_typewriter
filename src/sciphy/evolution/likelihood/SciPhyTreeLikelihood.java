package sciphy.evolution.likelihood;


import java.util.*;

import beast.base.core.Description;
import beast.base.core.Log;
import beast.base.evolution.likelihood.GenericTreeLikelihood;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.State;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.branchratemodel.BranchRateModel;
import beast.base.evolution.branchratemodel.StrictClockModel;
import beast.base.evolution.sitemodel.SiteModel;
import beast.base.evolution.substitutionmodel.SubstitutionModel;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeInterface;
import sciphy.evolution.substitutionmodel.SciPhySubstitutionModel;
import static sciphy.util.LogSum.logSum;

@Description("tree likelihood for a SciPhy alignment given a generic SiteModel, " +
        "a beast tree and a branch rate model. This is a version of the SciPhy likelihood using caching without a likelihoodCore implementation ")

public class SciPhyTreeLikelihood extends GenericTreeLikelihood {

    private static List<Integer> missingState;
    final public Input<RealParameter> originTimeInput = new Input<>("origin", "Duration of the experiment");

    final public Input<IntegerParameter> arrayLengthInput = new Input<>("arrayLength", "Number of positions in the target BC", Validate.REQUIRED);

    final public Input<Boolean> useScalingInput = new Input<Boolean>("useScaling", "Whether or not to scale the log likelihood", false,
            Validate.OPTIONAL);

    protected SciPhySubstitutionModel substitutionModel;
    protected BranchRateModel.Base branchRateModel;
    protected SiteModel.Base m_siteModel;
    protected double originTime;
    protected int nodeCount;
    protected int arrayLength;


    /**
     * flag to indicate the
     * // when CLEAN=0, nothing needs to be recalculated for the node
     * // when DIRTY=1 indicates a node partial needs to be recalculated
     * // when FILTHY=2 indicates the indices for the node need to be recalculated
     * // (often not necessary while node partial recalculation is required)
     */
    protected int hasDirt;

    /**
     * Lengths of the branches in the tree associated with each of the nodes
     * in the tree through their node  numbers. By comparing whether the
     * current branch length differs from stored branch lengths, it is tested
     * whether a node is dirty and needs to be recomputed (there may be other
     * reasons as well...).
     * These lengths take branch rate models in account.
     */
    protected double[] m_branchLengths;
    protected double[] storedBranchLengths;

    //to be able to have current/stored states in an analog way to the partials array, ancestral states are accessed/added
    //states with key being the hashcode of: nodeIndex + currentStatesIndex[nodeIndex] + nodeIndex
    public Hashtable<Integer, List<List<Integer>>> ancestralStates;
    public double[][][] partialLikelihoods;
    public double[] categoryLogLikelihoods;
    protected double[][] scalingFactors;
    protected boolean useScaling = false;


    private double scalingThreshold = 1.0E-100;

    protected int[] currentPartialsIndex;
    protected int[] storedPartialsIndex;

    protected int[] currentStatesIndex;
    protected int[] storedStatesIndex;


    @Override
    public void initAndValidate() {
        arrayLength = arrayLengthInput.get().getValue();
        if (arrayLength < 1 || (dataInput.get().getSiteCount() != arrayLength)) {
            throw new IllegalArgumentException(String.format(
                    "Invalid array length: Ensure that length >= 1 and matches alignment "));
        }
        nodeCount = treeInput.get().getNodeCount();
        if (nodeCount <= 2) {
            throw new IllegalArgumentException(String.format(
                    "Invalid tree input: single node/branch. Ensure that #nodes>2 "));
        }
        m_siteModel = (SiteModel.Base) siteModelInput.get();
        categoryLogLikelihoods = new double[m_siteModel.getCategoryCount()];
        m_siteModel.setDataType(dataInput.get().getDataType());
        substitutionModel = (SciPhySubstitutionModel) m_siteModel.substModelInput.get();

        m_branchLengths = new double[nodeCount];
        storedBranchLengths = new double[nodeCount];

        //TODO check that state count from alignment (i.e. data type) and substitution model are the same
        ancestralStates = new Hashtable<>();
        partialLikelihoods = new double[2][nodeCount][];

        currentPartialsIndex = new int[nodeCount];
        storedPartialsIndex = new int[nodeCount];

        currentStatesIndex = new int[nodeCount];
        storedStatesIndex = new int[nodeCount];

        if (branchRateModelInput.get() != null) {
            branchRateModel = branchRateModelInput.get();
        } else {
            branchRateModel = new StrictClockModel();
        }
        originTime = 0.0;
        if (originTimeInput.get() != null) {
            originTime = originTimeInput.get().getValue();
            if (originTime < 0.0) {
                throw new IllegalArgumentException(String.format(
                        "Invalid origin time input: ensure that origin>0"));
            }
        }

        if (useScalingInput.get()) {
            useScaling = true;
            scalingFactors = new double[2][nodeCount];
        }


        hasDirt = Tree.IS_FILTHY;

        // creating the missing state, encoded as an array of -1.
        missingState = new ArrayList<Integer>(){{
            add(-1);
            add(-1);
            add(-1);
            add(-1);
            add(-1);
        }};

        for (int i = 0; i < treeInput.get().getLeafNodeCount(); i++) {
            initLeafAncestors(i);
        }

        for (int i = 0; i < treeInput.get().getLeafNodeCount(); i++) {
            initLeafPartials(i);
        }

    }

    public SubstitutionModel getSubstitutionModel() {
        return substitutionModel;
    }

    @Override
    public List<String> getArguments() {
        return null;
    }

    @Override
    public List<String> getConditions() {
        return null;
    }

    @Override
    public void sample(State state, Random random) {
    }

    @Override
    public double calculateLogP() {
        final TreeInterface tree = treeInput.get();

        if(originTime != 0.0) {
            if (tree.getRoot().getHeight() >= originTime) {
                return Double.NEGATIVE_INFINITY;
            }
        }

        for (int i = 0; i < m_siteModel.getCategoryCount(); i++) {
            //adjust clock rate for the given category
            traverse(tree.getRoot(), i);

            if (originTime == 0.0) {
                //sum of all partial likelihoods at the root
                int rootNr = tree.getRoot().getNr();
                categoryLogLikelihoods[i] = Math.log(Arrays.stream(partialLikelihoods[currentPartialsIndex[rootNr]][rootNr]).sum()) + getLogScalingFactor();
            } else {
                //the tree log likelihood is the log(p) of unedited state at the origin
                categoryLogLikelihoods[i] = Math.log(calculateOriginPartial(tree.getRoot(), i)) + getLogScalingFactor();

            }
        }
        logP = logSum(categoryLogLikelihoods, categoryLogLikelihoods.length) - Math.log(m_siteModel.getCategoryCount());
        return logP;
    }

    /**
     * Scale the partials at a given node. This uses a scaling suggested by Ziheng Yang in
     * Yang (2000) J. Mol. Evol. 51: 423-432
     * <p/>
     * This function looks over the partial likelihoods for each state at each pattern
     * and finds the largest. If this is less than the scalingThreshold (currently set
     * to 1E-40) then it rescales the partials for that pattern by dividing by this number
     * (i.e., normalizing to between 0, 1). It then stores the log of this scaling.
     * This is called for every internal node after the partials are calculated so provides
     * most of the performance hit. Ziheng suggests only doing this on a proportion of nodes
     * but this sounded like a headache to organize (and he doesn't use the threshold idea
     * which improves the performance quite a bit).
     *
     * @param nodeNumber
     */
    protected void scalePartials(int nodeNumber) {

        double scaleFactor = 0.0;

        //find the highest partial likelihood
        //is node number same as nodeIndex
        for (int k = 0; k < partialLikelihoods[currentPartialsIndex[nodeNumber]][nodeNumber].length; k++) {

            if (partialLikelihoods[currentPartialsIndex[nodeNumber]][nodeNumber][k] > scaleFactor) {
                scaleFactor = partialLikelihoods[currentPartialsIndex[nodeNumber]][nodeNumber][k];
            }

        }
        //if this partial is smaller than the threshold, scale the partials
        if (scaleFactor < scalingThreshold) {

            for (int k = 0; k < partialLikelihoods[currentPartialsIndex[nodeNumber]][nodeNumber].length; k++) {
                partialLikelihoods[currentPartialsIndex[nodeNumber]][nodeNumber][k] /= scaleFactor;
            }
            // save the log(scaling factors)
            scalingFactors[currentPartialsIndex[nodeNumber]][nodeNumber] = Math.log(scaleFactor);

        } else {
            scalingFactors[currentPartialsIndex[nodeNumber]][nodeNumber] = 0.0;
        }

    }

    /**
     * Calculate partial likelihoods for a given leaf node, and fill the corresponding partialLikelihood array
     */
    protected void initLeafPartials(int nodeNr) {

        double[] leafPartialLikelihoods = initPartialLikelihoodsLeaf(ancestralStates.get(makeCachingIndexStates(nodeNr)).size());
        this.partialLikelihoods[0][nodeNr] = new double[leafPartialLikelihoods.length];
        this.partialLikelihoods[1][nodeNr] = new double[leafPartialLikelihoods.length];
        System.arraycopy(leafPartialLikelihoods, 0, this.partialLikelihoods[0][nodeNr], 0, leafPartialLikelihoods.length);

    }


    public int makeCachingIndexStates(int nodeIndex) {
        int node = nodeIndex + 1;
        String forHashing = node + "" +  currentStatesIndex[nodeIndex] + ""+ node;
        return forHashing.hashCode();

    }


    /**
     * Calculate the set of ancestral states for a given leaf node, and fill the corresponding AncestralStates hashmap
     */
    protected void initLeafAncestors(int nodeNr) {

        List<List<Integer>> possibleLeafAncestors = getPossibleAncestors(dataInput.get().getCounts().get(nodeNr));
        ancestralStates.put(makeCachingIndexStates(nodeNr), possibleLeafAncestors);

    }


    /**
     * This implements a postorder traversal of the tree to fill the ancestralStates hashmap and corresponding partialLikelihood array.
     */
    protected int traverse(Node node, int categoryId) {

        int update = (node.isDirty() | hasDirt);
        int nodeIndex = node.getNr();

        final double branchRate = branchRateModel.getRateForBranch(node);
        final double branchTime = node.getLength() * branchRate;

        if (!node.isRoot() && (update != Tree.IS_CLEAN || branchTime != m_branchLengths[nodeIndex])) {
            m_branchLengths[nodeIndex] = branchTime;
            update |= Tree.IS_DIRTY;
        }

        if (!node.isLeaf()) {

            final Node child1 = node.getLeft();
            final int update1 = traverse(child1, categoryId);
            final Node child2 = node.getRight();
            final int update2 = traverse(child2, categoryId);

            // If either child node was updated then update this node too
            if (update1 != Tree.IS_CLEAN || update2 != Tree.IS_CLEAN) {

                update |= (update1 | update2);

                if (update >= Tree.IS_FILTHY) {
                    setNodeStatesForUpdate(nodeIndex);
                    calculateStates(nodeIndex, child1.getNr(), child2.getNr());
                }

                setNodePartialsForUpdate(nodeIndex);
                calculatePartials(nodeIndex, child1, child2, categoryId);

                if (useScaling) {
                    scalePartials(nodeIndex);
                }

            }
        }

        return update;
    }

    /**
     * Construct a set of possible ancestral states at an internal node by intersection of children sets, updates the
     * AncestralStates hashmap with the resulting set.
     */
    public void calculateStates(int nodeNr, int child1Nr, int child2Nr) {

        List<List<Integer>> ancSetChild1 = ancestralStates.get(makeCachingIndexStates(child1Nr));
        List<List<Integer>> ancSetChild2 = ancestralStates.get(makeCachingIndexStates(child2Nr));



        // create the set of ancestral states for the missing state
        List<List<Integer>> ancestorsMissing = new ArrayList();
        ancestorsMissing.add(missingState);

        // both sets are missing states, the intersection is the missing state
        //todo rearrage these statements
        if( ancSetChild2.equals(ancestorsMissing) && ancSetChild1.equals(ancestorsMissing) ) {
            ancestralStates.put(makeCachingIndexStates(nodeNr), ancestorsMissing);
        }
        else {


            if (ancSetChild1.equals(ancestorsMissing)) {
                //the intersection is the child2 set
                List<List<Integer>> ancSetNode = new ArrayList<>(ancSetChild2);
                ancestralStates.put(makeCachingIndexStates(nodeNr), ancSetNode);
            }
            else if (ancSetChild2.equals(ancestorsMissing)) {
                //the intersection is the child1 set
                List<List<Integer>> ancSetNode = new ArrayList<>(ancSetChild1);
                ancestralStates.put(makeCachingIndexStates(nodeNr), ancSetNode);
            }
            else {
                //intersection needs to be computed
                List<List<Integer>> ancSetNode = new ArrayList<>(ancSetChild1);
                ancSetNode.retainAll(ancSetChild2);
                ancestralStates.put(makeCachingIndexStates(nodeNr), ancSetNode);
            }

        }


    }

    public void setNodePartialsForUpdate(int nodeIndex) {
        currentPartialsIndex[nodeIndex] = 1 - currentPartialsIndex[nodeIndex];
    }

    public void setNodeStatesForUpdate(int nodeIndex) {
        currentStatesIndex[nodeIndex] = 1 - currentStatesIndex[nodeIndex];
    }


    /**
     * This function calculates partial likelihoods for all possible states at a node given its children partials
     * and sets the corresponding partial likelihoods, for all possible states at node nodeNr
     */
    public void calculatePartials(int nodeNr, Node child1, Node child2, int categoryId) {

        //initialize an array for the partials
        double[] partials = new double[ancestralStates.get(makeCachingIndexStates(nodeNr)).size()];

        for (int stateIndex = 0; stateIndex < ancestralStates.get(makeCachingIndexStates(nodeNr)).size(); ++stateIndex) {

            List<Integer> startState = ancestralStates.get(makeCachingIndexStates(nodeNr)).get(stateIndex);

            double child1PartialLikelihoodState = calculatePartialLikelihoodState(startState, child1, categoryId);
            double child2PartialLikelihoodState = calculatePartialLikelihoodState(startState, child2, categoryId);

            partials[stateIndex] = child1PartialLikelihoodState * child2PartialLikelihoodState;
        }

        partialLikelihoods[currentPartialsIndex[nodeNr]][nodeNr] = partials;

    }

    /**
     * This function calculates the likelihood of the unedited state at the origin given partial likelihoods at the root
     * node
     *
     * @return likelihood of the unedited barcode at t = origin
     */

    public double calculateOriginPartial(Node rootNode, int categoryId) {

        //the start state is the unedited sciphy barcode
        List<Integer> startState = Arrays.asList(0, 0, 0, 0, 0);
        double partialAtOrigin = calculatePartialLikelihoodState(startState, rootNode, categoryId);
        return partialAtOrigin;

    }

    /**
     * This function calculates the partial likelihood term of a specific state at a node derived on a branch leading to
     * a child node
     *
     * @return partial likelihood for a state at a node given partials at a node childNode
     */
    public double calculatePartialLikelihoodState(List<Integer> startState, Node childNode, int categoryId) {

        final double branchRate = branchRateModel.getRateForBranch(childNode);
        double statePartialLikelihood = 0;
        final double jointBranchRate = m_siteModel.getRateForCategory(categoryId, childNode) * branchRate;
        double distance;

        //initialise evolutionary distance
        if (childNode.isRoot()) {
            distance = (originTime - childNode.getHeight()) * jointBranchRate;
        } else {
            distance = childNode.getLength() * jointBranchRate;
        }
        // calculate partials
        if (childNode.isLeaf()) {

            List<Integer> endState = ancestralStates.get(makeCachingIndexStates(childNode.getNr())).get(0);
            statePartialLikelihood += substitutionModel.getSequenceTransitionProbability(startState, endState, distance, this.arrayLength);

        } else {

            for (int endStateIndex = 0; endStateIndex < ancestralStates.get(makeCachingIndexStates(childNode.getNr())).size(); ++endStateIndex) {

                List<Integer> endState = ancestralStates.get(makeCachingIndexStates(childNode.getNr())).get(endStateIndex);

                // if the end state has non-null partial likelihood
                if (partialLikelihoods[currentPartialsIndex[childNode.getNr()]][childNode.getNr()][endStateIndex] != 0.0) {

                    statePartialLikelihood = statePartialLikelihood + substitutionModel.getSequenceTransitionProbability(startState, endState, distance, this.arrayLength) *
                            partialLikelihoods[currentPartialsIndex[childNode.getNr()]][childNode.getNr()][endStateIndex];

                }
            }
        }
        return statePartialLikelihood;
    }

    /**
     * This function initialises an array of partial likelihoods for a leaf node, the partial likelihood is 1 for
     * the observed sequence and 0 for everything else. The size corresponds to the total number of possible ancestral states.
     *
     * @return array of partial likelihoods at leaf node
     */
    public double[] initPartialLikelihoodsLeaf(int size) {

        double[] leafPartials = new double[size];
        leafPartials[0] = 1;
        return leafPartials;
    }

    /**
     * This function returns all possible ancestral states given a sequence.
     * Because sciphy sequences record ordered edits, ancestral states are obtained by sequentially removing edits
     * along the sequence (from any insert (1 to N)  to 0)
     *
     * @return a list of possible ancestral sciphy barcode states
     */
    public static List<List<Integer>> getPossibleAncestors(List<Integer> sequence) {

        List<List<Integer>> ancestors = new ArrayList();
        ancestors.add(sequence);

        //possible ancestors of a missing state can be any state, we return the missing state itself as a Wild Card
        if(sequence.equals(missingState)) {
            return ancestors;
        }

        else {
            List<Integer> ancestor = new ArrayList<>(sequence);
            for (int i = sequence.size() - 1; i >= 0; --i) {
                if (sequence.get(i) != 0) {
                    ancestor.set(i, 0);
                    ancestors.add(new ArrayList<>(ancestor));
                }
            }
            return ancestors;
        }
    }

    /**
     * This function returns the scaling factor for that pattern by summing over
     * the log scalings used at each node. If scaling is off then this just returns
     * a 0.
     *
     * @return the log scaling factor
     */
    public double getLogScalingFactor() {

        double logScalingFactor = 0.0;
        if (useScaling) {
            for (int i = 0; i < nodeCount; i++) {
                logScalingFactor += scalingFactors[currentPartialsIndex[i]][i];
            }
        }
        return logScalingFactor;
    }


    /**
     * check state for changed variables and update temp results if necessary *
     */
    //requires recalculation if the data has changed, if the sitemodel (or sub model has changed), the branch rate model, or the tree has changed.
    @Override
    protected boolean requiresRecalculation() {
        hasDirt = Tree.IS_CLEAN;

        if (dataInput.get().isDirtyCalculation()) {
            hasDirt = Tree.IS_FILTHY;
            return true;
        }
        if (m_siteModel.isDirtyCalculation()) {
            hasDirt = Tree.IS_DIRTY;
            return true;
        }
        if (branchRateModel != null && branchRateModel.isDirtyCalculation()) {
            //m_nHasDirt = Tree.IS_DIRTY;
            return true;
        }
        return treeInput.get().somethingIsDirty();
    }

    @Override
    public void store() {

        super.store();
        System.arraycopy(m_branchLengths, 0, storedBranchLengths, 0, m_branchLengths.length);
        System.arraycopy(currentPartialsIndex, 0, storedPartialsIndex, 0, nodeCount);
        System.arraycopy(currentStatesIndex, 0, storedStatesIndex, 0, nodeCount);
    }

    //TODO do we need unstore??? We think we don't because when scaling is active, it is for the entire likelihood

    @Override
    public void restore() {

        super.restore();
        double[] tmp = m_branchLengths;
        m_branchLengths = storedBranchLengths;
        storedBranchLengths = tmp;

        int[] tmp2 = currentPartialsIndex;
        currentPartialsIndex = storedPartialsIndex;
        storedPartialsIndex = tmp2;

        int[] tmp3 = currentStatesIndex;
        currentStatesIndex = storedStatesIndex;
        storedStatesIndex = tmp3;
    }


}
