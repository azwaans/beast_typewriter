package lineageTree.distributions;


import java.util.*;

import beast.core.Description;
import beast.core.Distribution;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.State;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Alignment;
import beast.evolution.branchratemodel.BranchRateModel;
import beast.evolution.branchratemodel.StrictClockModel;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.sitemodel.SiteModelInterface;
import beast.evolution.substitutionmodel.SubstitutionModel;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.evolution.tree.TreeInterface;
import lineageTree.substitutionmodel.TypewriterSubstitutionModel;

import static java.lang.Math.log1p;

@Description("tree likelihood for a Typewriter alignment given a generic SiteModel, " +
        "a beast tree and a branch rate model. This is a version of the Typewriter likelihood using caching without a likelihoodCore implementation ")

public class TypewriterTreeLikelihood extends Distribution {

    final public Input<Alignment> dataInput = new Input<>("data", "sequence data for the beast.tree", Validate.REQUIRED);

    final public Input<TreeInterface> treeInput = new Input<>("tree", "phylogenetic beast.tree with sequence data in the leafs", Validate.REQUIRED);

    final public Input<SiteModelInterface> siteModelInput = new Input<>("siteModel", "site model for leafs in the beast.tree", Validate.REQUIRED);

    final public Input<BranchRateModel.Base> branchRateModelInput = new Input<>("branchRateModel",
            "A model describing the rates on the branches of the beast.tree.");

    final public Input<RealParameter> originTimeInput = new Input<>("origin", "Duration of the experiment");

    final public Input<IntegerParameter> arrayLengthInput = new Input<>("arrayLength", "Number of positions in the target BC");

    final public Input<Boolean> useScalingInput = new Input<Boolean>("useScaling", "Whether or not to scale the log likelihood", false,
            Validate.OPTIONAL);

    protected TypewriterSubstitutionModel substitutionModel;
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

    //to be able to have current/stored states in an analog way to the partials array, ancestral states are we get/add
    //states with key : (NodeNr + 1) + (current ? 0:1) * (NodeNr+1)
    public Hashtable<Integer,List<List<Integer>>> ancestralStates ;
    public double[][][] partialLikelihoods ;
    public double[] categoryLogLikelihoods ;
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
        if(arrayLength < 1 || (dataInput.get().getSiteCount() != arrayLength)) {
            throw new IllegalArgumentException(String.format(
                    "Invalid array length: Ensure that length >= 1 and matches alignment "));


        }
        nodeCount = treeInput.get().getNodeCount();
        if(nodeCount <=2) {
            throw new IllegalArgumentException(String.format(
                    "Invalid tree input: single node/branch. Ensure that #nodes>2 "));
        }
        m_siteModel = (SiteModel.Base) siteModelInput.get();
        categoryLogLikelihoods = new double[m_siteModel.getCategoryCount()];
        m_siteModel.setDataType(dataInput.get().getDataType());

        substitutionModel = (TypewriterSubstitutionModel)  m_siteModel.substModelInput.get();
        substitutionModel.setTargetBClength(arrayLength);

        m_branchLengths = new double[nodeCount];
        storedBranchLengths = new double[nodeCount];


        //TODO check that state count from alignment (i.e. data type) and substitution model are the same
        //TODO INITIALISE EVERYTHING TO THE NUMBER OF NODES + MAX NUMBER OF STATES.
        ancestralStates = new Hashtable<>() ;
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

        if(useScalingInput.get()){
            useScaling = true;
            scalingFactors = new double[2][nodeCount];
        }

        //initialize has dirt flag (this is updated by requires recalcuations):
        hasDirt = Tree.IS_FILTHY;

        for (int i=0; i< treeInput.get().getLeafNodeCount(); i++) {
            initLeafAncestors(i);
        }

        for (int i=0; i< treeInput.get().getLeafNodeCount(); i++) {
            initLeafPartials(i);
        }

    }

    public SubstitutionModel getSubstitutionModel() {return substitutionModel;}

    @Override
    public List<String> getArguments() {return null;}

    @Override
    public List<String> getConditions() {return null;}

    @Override
    public void sample(State state, Random random) {}

    @Override
    public double calculateLogP() {

        final TreeInterface tree = treeInput.get();
        substitutionModel.setTargetBClength(arrayLength);

        for (int i = 0; i < m_siteModel.getCategoryCount(); i++) {
            //adjust clock rate for the given category
            traverseLikelihood(tree.getRoot(),i);

            if (originTime == 0.0) {
                //sum of all partial likelihoods at the root
                int rootNr = tree.getRoot().getNr();
                categoryLogLikelihoods[i] = Math.log(Arrays.stream(partialLikelihoods[currentPartialsIndex[rootNr]][rootNr]).sum()) + getLogScalingFactor();
            } else {
                //the tree log likelihood is the log(p) of unedited state at the origin
                categoryLogLikelihoods[i] = Math.log(calculateOriginPartial(tree.getRoot(),i)) + getLogScalingFactor();

            }
        }
        logP =  log_sum(categoryLogLikelihoods, categoryLogLikelihoods.length) - Math.log(m_siteModel.getCategoryCount());
        return logP;
    }

    double log_sum(double la[], int num_elements)
    {
        // Assume index_of_max() finds the maximum element
        // in the array and returns its index
        double max = la[0];
        int index = 0;
        for (int i = 0; i < la.length; i++)
        {
            if (max < la[i])
            {
                max = la[i];
                index = i;
            }
        }

        double sum_exp = 0;
        for (int i = 0; i < num_elements; i++) {
            if (i == index) {
                continue;
            }
            sum_exp += Math.exp(la[i] - la[index]);
        }

        return la[index] + log1p(sum_exp);
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
    //TODO be careful: if consecutive scaling was added: needed store and restore!
    protected void scalePartials(int nodeNumber) {

        double scaleFactor = 0.0;

        //find the highest partial likelihood
        //is node number same as nodeIndex
        for (int k = 0; k < partialLikelihoods[currentPartialsIndex[nodeNumber]][nodeNumber].length; k++) {

            if(partialLikelihoods[currentPartialsIndex[nodeNumber]][nodeNumber][k] > scaleFactor) {
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
     * This function implements a postorder traversal of the tree to obtain all possible sets of ancestral state sets. This is only used for testing purposes
     */

    public void traverseAncestral(Node node) {

       if(!node.isLeaf()) {

            final Node child1 = node.getLeft();
            final Node child2 = node.getRight();

            traverseAncestral(child1);
            traverseAncestral(child2);

            List<List<Integer>> ancSetChild1 = ancestralStates.get(child1.getNr() );
            List<List<Integer>> ancSetChild2 = ancestralStates.get(child2.getNr() );

           // intersection of children ancestral states
            List<List<Integer>> ancSetNode = new ArrayList<>(ancSetChild1);
            ancSetNode.retainAll(ancSetChild2);

            ancestralStates.put(node.getNr()+1, ancSetNode);

       }

    }



    protected void initLeafPartials(int nodeNr) {

        //create and fill tip partials
        double[] leafPartialLikelihoods = initPartialLikelihoodsLeaf(ancestralStates.get((nodeNr+1) + currentStatesIndex[nodeNr]*(nodeNr +1)).size());
        this.partialLikelihoods[0][nodeNr] = new double[leafPartialLikelihoods.length];
        this.partialLikelihoods[1][nodeNr] = new double[leafPartialLikelihoods.length];
        System.arraycopy(leafPartialLikelihoods, 0, this.partialLikelihoods[0][nodeNr], 0, leafPartialLikelihoods.length);

    }

    protected void initLeafAncestors(int nodeNr) {

        List<List<Integer>> possibleLeafAncestors = getPossibleAncestors(dataInput.get().getCounts().get(nodeNr));
        ancestralStates.put(nodeNr+1, possibleLeafAncestors);

    }


    /**
     * This function implements a postorder traversal of the tree to fill the partialLikelihood array
     *
     */
    protected int traverseLikelihood(Node node, int categoryId) {

        int update = (node.isDirty() | hasDirt);
        int nodeIndex = node.getNr();

        final double branchRate = branchRateModel.getRateForBranch(node);
        final double branchTime = node.getLength() * branchRate;

        if (!node.isRoot() && (update != Tree.IS_CLEAN || branchTime != m_branchLengths[nodeIndex])) {
            m_branchLengths[nodeIndex] = branchTime;
            update |= Tree.IS_DIRTY;
        }

        if(!node.isLeaf()) {

            final Node child1 = node.getLeft();
            final int update1 = traverseLikelihood(child1, categoryId);
            final Node child2 = node.getRight();
            final int update2 = traverseLikelihood(child2, categoryId);

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


    public void calculateStates(int nodeNr, int child1Nr, int child2Nr) {

        List<List<Integer>> ancSetChild1 = ancestralStates.get(child1Nr +1 + currentStatesIndex[child1Nr]*(child1Nr+1));
        List<List<Integer>> ancSetChild2 = ancestralStates.get(child2Nr +1 + currentStatesIndex[child2Nr]*(child2Nr+1));

        List<List<Integer>> ancSetNode = new ArrayList<>(ancSetChild1);
        // intersection of children ancestral states
        ancSetNode.retainAll(ancSetChild2);

        ancestralStates.put((nodeNr +1) + currentStatesIndex[nodeNr]*(nodeNr+1), ancSetNode);
    }

    public void setNodePartialsForUpdate(int nodeIndex) {
        currentPartialsIndex[nodeIndex] = 1 - currentPartialsIndex[nodeIndex];
    }

    public void setNodeStatesForUpdate(int nodeIndex) {
        currentStatesIndex[nodeIndex] = 1 - currentStatesIndex[nodeIndex];
    }



    /**
     * This function calculates partial likelihoods for all possible states at node given its children partials
     *
     * @return partial likelihoods for states at node nodeNr
     */
    public void calculatePartials(int nodeNr, Node child1, Node child2, int categoryId ) {

        //initialize an array for the partials
        double[] partials = new double[ancestralStates.get((nodeNr+1) + currentStatesIndex[nodeNr]*(nodeNr +1)).size()];

            for (int stateIndex = 0; stateIndex < ancestralStates.get((nodeNr +1) + currentStatesIndex[nodeNr]*(nodeNr+1)).size(); ++stateIndex) {
                
                List<Integer> startState = ancestralStates.get((nodeNr+1) + currentStatesIndex[nodeNr]*(nodeNr+1)).get(stateIndex);
                
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

        //the start state is the unedited typewriter barcode
        List<Integer> startState = Arrays.asList(0,0,0,0,0);
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

            List<Integer> endState = ancestralStates.get(childNode.getNr()+1 + currentStatesIndex[childNode.getNr()]*(childNode.getNr()+1)).get(0);
            statePartialLikelihood += substitutionModel.getSequenceTransitionProbability(startState, endState, distance);

        } else {

            for (int endStateIndex = 0; endStateIndex < ancestralStates.get(childNode.getNr()+1 + currentStatesIndex[childNode.getNr()]*(childNode.getNr()+1)).size(); ++endStateIndex) {

                List<Integer> endState = ancestralStates.get(childNode.getNr() +1  + currentStatesIndex[childNode.getNr()]*(childNode.getNr()+1)).get(endStateIndex);

                // if the end state has non null partial likelihood
                if (partialLikelihoods[currentPartialsIndex[childNode.getNr()]][childNode.getNr()][endStateIndex] != 0.0) {

                    statePartialLikelihood = statePartialLikelihood + substitutionModel.getSequenceTransitionProbability(startState, endState, distance) *
                            partialLikelihoods[currentPartialsIndex[childNode.getNr()]][childNode.getNr()][endStateIndex];

                }
            }
        }
        return statePartialLikelihood ;
    }

    /**
     * This function initialises an array of partial likelihoods for a leaf node, the partial likelihood is 1 for
     * the observed sequence and 0 for everything else. The size corresponds to the total number of possible ancestral states.
     *
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
     * Because typewriter sequences record ordered edits, ancestral states are obtained by sequentially removing edits
     * along the sequence (from any insert (1 to N)  to 0)
     *
     * @return a list of possible ancestral typewriter barcode states
     */
    public static List<List<Integer>> getPossibleAncestors(List<Integer> sequence) {

        List<List<Integer>> ancestors = new ArrayList();
        ancestors.add(sequence);

        List<Integer> ancestor = new ArrayList<>(sequence);
        for(int i = sequence.size()-1; i >= 0; --i) {
            if(sequence.get(i) != 0) {
                ancestor.set(i, 0);
                ancestors.add(new ArrayList<>(ancestor));
            }
        }
        return ancestors;
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
//        Log.info.println("store");
//        if (likelihoodCore != null) {
//            likelihoodCore.store();
//        }
        super.store();
        System.arraycopy(m_branchLengths, 0, storedBranchLengths, 0, m_branchLengths.length);
        System.arraycopy(currentPartialsIndex, 0, storedPartialsIndex, 0, nodeCount);
        System.arraycopy(currentStatesIndex, 0, storedStatesIndex, 0, nodeCount);
    }

    //TODO do we need unstore??? 
//    public void unstore() {
//        System.arraycopy(storedPartialsIndex, 0, currentPartialsIndex, 0, nodeCount);
//    }

    @Override
    public void restore() {

//        Log.info.println("restore");

//        if (likelihoodCore != null) {
//            likelihoodCore.restore();
//        }
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
