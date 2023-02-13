package lineageTree.distributions;


import java.util.*;

import beast.core.Description;
import beast.core.Distribution;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.State;
import beast.core.parameter.RealParameter;
import beast.core.util.Log;
import beast.evolution.alignment.Alignment;
import beast.evolution.branchratemodel.BranchRateModel;
import beast.evolution.branchratemodel.StrictClockModel;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.sitemodel.SiteModelInterface;
import beast.evolution.substitutionmodel.SubstitutionModel;
import beast.evolution.tree.Node;
import beast.evolution.tree.TreeInterface;
import lineageTree.substitutionmodel.TypewriterSubstitutionModelHomogeneous;

import static java.lang.Math.log1p;


@Description("tree likelihood for a Typewriter alignment given a generic SiteModel, " +
        "a beast tree and a branch rate model. This is a brute-force approach with no use of BEAST treelikelihood architecture ")

public class TypewriterTreeLikelihood extends Distribution {

    final public Input<Alignment> dataInput = new Input<>("data", "sequence data for the beast.tree", Validate.REQUIRED);

    final public Input<TreeInterface> treeInput = new Input<>("tree", "phylogenetic beast.tree with sequence data in the leafs", Validate.REQUIRED);

    final public Input<SiteModelInterface> siteModelInput = new Input<>("siteModel", "site model for leafs in the beast.tree", Validate.REQUIRED);

    final public Input<BranchRateModel.Base> branchRateModelInput = new Input<>("branchRateModel",
            "A model describing the rates on the branches of the beast.tree.");

    final public Input<RealParameter> originTimeInput = new Input<>("origin", "Duration of the experiment");

    final public Input<Boolean> useScalingInput = new Input<Boolean>("useScaling", "Whether or not to scale the log likelihood", false,
            Validate.OPTIONAL);

    protected TypewriterSubstitutionModelHomogeneous substitutionModel;
    protected BranchRateModel.Base branchRateModel;
    protected SiteModel.Base m_siteModel;
    protected double[] m_branchLengths;
    protected double originTime;
    protected int nodeCount;


    public Hashtable<Integer,List<List<Integer>>> ancestralStates ;
    public double[][] partialLikelihoods ;
    public double[] categoryLogLikelihoods ;
    protected double[] scalingFactors;
    protected boolean useScaling = false;


    private double scalingThreshold = 1.0E-100;



    @Override
    public void initAndValidate() {

        nodeCount = treeInput.get().getNodeCount();
        m_siteModel = (SiteModel.Base) siteModelInput.get();
        categoryLogLikelihoods = new double[m_siteModel.getCategoryCount()];
        m_siteModel.setDataType(dataInput.get().getDataType());

        substitutionModel = (TypewriterSubstitutionModelHomogeneous)  m_siteModel.substModelInput.get();
        m_branchLengths = new double[nodeCount];
        ancestralStates = new Hashtable<>() ;

        //TODO check that state count from alignment (i.e. data type) and substitution model are the same

        if (branchRateModelInput.get() != null) {
            branchRateModel = branchRateModelInput.get();
        } else {
            branchRateModel = new StrictClockModel();
        }
        originTime = 0.0;
        if (originTimeInput.get() != null) {
            originTime = originTimeInput.get().getValue();
        }

        partialLikelihoods = new double[nodeCount][];



        if(useScalingInput.get()){
            useScaling = true;
            scalingFactors = new double[nodeCount];
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

        //1st step : calculate all ancestral states in a Postorder traversal
        //TODO potentially find other data structure than hashmaps
        final TreeInterface tree = treeInput.get();
        traverseAncestral(tree.getRoot());

        //2nd step : calculate likelihood with these ancestral states
        // size of the partial likelihoods at each node = nr of the ancestral states at that node.


        for (int i = 0; i < m_siteModel.getCategoryCount(); i++) {
            //adjust clock rate for the given category
            traverseLikelihood(tree.getRoot(),i);

            if (originTime == 0.0) {
                //sum of all partial likelihoods at the root
                categoryLogLikelihoods[i] = Math.log(Arrays.stream(partialLikelihoods[tree.getRoot().getNr()]).sum()) + getLogScalingFactor();
            } else {
                //the tree log likelihood is the log(p) of unedited state at the origin
                //TODO check that origin logP is calculated
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
    protected void scalePartials(int nodeNumber) {

        double scaleFactor = 0.0;

        //find the highest partial likelihood
        //is node number same as nodeIndex
        for (int k = 0; k < partialLikelihoods[nodeNumber].length; k++) {
            if(partialLikelihoods[nodeNumber][k] > scaleFactor) {
                scaleFactor = partialLikelihoods[nodeNumber][k];
            }

        }
        //if this partial is smaller than the threshold, scale the partials
        if (scaleFactor < scalingThreshold) {

            for (int k = 0; k < partialLikelihoods[nodeNumber].length; k++) {

                partialLikelihoods[nodeNumber][k] /= scaleFactor;

            }
            // save the log(scaling factors)
            scalingFactors[nodeNumber] = Math.log(scaleFactor);

        } else {
            scalingFactors[nodeNumber] = 0.0;
        }

    }

    /**
     * This function implements a postorder traversal of the tree to obtain all possible sets of ancestral state sets
     */
    public void traverseAncestral(Node node) {

        if (node.isLeaf() ) {
            List<List<Integer>> possibleLeafAncestors = getPossibleAncestors(dataInput.get().getCounts().get(node.getNr()));
            ancestralStates.put(node.getNr(), possibleLeafAncestors);

        } else {

            final Node child1 = node.getLeft();
            final Node child2 = node.getRight();

            traverseAncestral(child1);
            traverseAncestral(child2);

            List<List<Integer>> ancSetChild1 = ancestralStates.get(child1.getNr());
            List<List<Integer>> ancSetChild2 = ancestralStates.get(child2.getNr());

            List<List<Integer>> ancSetNode = new ArrayList<>(ancSetChild1);
            // intersection of children ancestral states
            ancSetNode.retainAll(ancSetChild2);

            ancestralStates.put(node.getNr(), ancSetNode);

        }


    }

    /**
     * This function implements a postorder traversal of the tree to fill the partialLikelihood array
     *
     */
    public void traverseLikelihood(Node node, int categoryId) {

        if( node != null ) {

            if (node.isLeaf()) {
                
                double[] leafPartialLikelihoods = initPartialLikelihoodsLeaf(ancestralStates.get(node.getNr()).size());
                partialLikelihoods[node.getNr()] = leafPartialLikelihoods;


            } else {

                final Node child1 = node.getLeft();
                final Node child2 = node.getRight();

                traverseLikelihood(child1, categoryId);
                traverseLikelihood(child2, categoryId);

                double[] partials = calculatePartials(node.getNr(),child1,child2,categoryId);
                partialLikelihoods[node.getNr()] = partials;

                if (useScaling) {
                    scalePartials(node.getNr());
                }

            }
        }



    }

    /**
     * This function calculates partial likelihoods for all possible states at node given its children partials
     *
     * @return partial likelihoods for states at node nodeNr
     */
    public double[] calculatePartials(int nodeNr, Node child1, Node child2, int categoryId ) {

        //initialize an array for the partials
        double[] partials = new double[ancestralStates.get(nodeNr).size()];

            for (int stateIndex = 0; stateIndex < ancestralStates.get(nodeNr).size(); ++stateIndex) {
                
                List<Integer> startState = ancestralStates.get(nodeNr).get(stateIndex);
                
                double child1PartialLikelihoodState = calculatePartialLikelihoodState(startState, child1, categoryId);
                double child2PartialLikelihoodState = calculatePartialLikelihoodState(startState, child2, categoryId);
                
                partials[stateIndex] = child1PartialLikelihoodState * child2PartialLikelihoodState;
            }
        return partials;

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
//        for (int i = 0; i < m_siteModel.getCategoryCount(); i++) {
            final double jointBranchRate = m_siteModel.getRateForCategory(categoryId, childNode) * branchRate;
           Log.info.println("site model rate at category"+"i = " + m_siteModel.getRateForCategory(categoryId, childNode));
            Log.info.println("joint branch rate" + jointBranchRate);


            double distance;

            //initialise evolutionary distance
            // TODO rename distance
            if (childNode.isRoot()) {
                distance = (originTime - childNode.getHeight()) * jointBranchRate;
            } else {
                distance = childNode.getLength() * jointBranchRate;
            }

            // calculate partials


            if (childNode.isLeaf()) {

                List<Integer> endState = ancestralStates.get(childNode.getNr()).get(0);
                statePartialLikelihood += substitutionModel.getSequenceTransitionProbability(startState, endState, distance);

            } else {

                for (int endStateIndex = 0; endStateIndex < ancestralStates.get(childNode.getNr()).size(); ++endStateIndex) {

                    List<Integer> endState = ancestralStates.get(childNode.getNr()).get(endStateIndex);

                    // if the end state has non null partial likelihood
                    if (partialLikelihoods[childNode.getNr()][endStateIndex] != 0.0) {

                        statePartialLikelihood = statePartialLikelihood + substitutionModel.getSequenceTransitionProbability(startState, endState, distance) *
                                partialLikelihoods[childNode.getNr()][endStateIndex];

                    }
                }
            }
//        }

//        return statePartialLikelihood / m_siteModel.getCategoryCount();
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
                logScalingFactor += scalingFactors[i];
            }
        }
        return logScalingFactor;
    }

}
