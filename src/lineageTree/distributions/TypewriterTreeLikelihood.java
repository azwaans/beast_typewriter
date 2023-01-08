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
import beast.evolution.tree.Tree;
import beast.evolution.tree.TreeInterface;
import lineageTree.substitutionmodel.TypewriterSubstitutionModel;
import lineageTree.substitutionmodel.TypewriterSubstitutionModelHomogeneous;


@Description("tree likelihood for a Typewriter alignment given a generic SiteModel, " +
        "a beast tree and a branch rate model. This is a brute-force approach with no use of BEAST treelikelihood architecture ")

public class TypewriterTreeLikelihood extends Distribution {

    final public Input<Alignment> dataInput = new Input<>("data", "sequence data for the beast.tree", Validate.REQUIRED);

    final public Input<TreeInterface> treeInput = new Input<>("tree", "phylogenetic beast.tree with sequence data in the leafs", Validate.REQUIRED);

    final public Input<SiteModelInterface> siteModelInput = new Input<>("siteModel", "site model for leafs in the beast.tree", Validate.REQUIRED);

    final public Input<BranchRateModel.Base> branchRateModelInput = new Input<>("branchRateModel",
            "A model describing the rates on the branches of the beast.tree.");

    final public Input<RealParameter> originTimeInput = new Input<>("origin", "Duration of the experiment");


    protected TypewriterSubstitutionModelHomogeneous substitutionModel;
    protected BranchRateModel.Base branchRateModel;
    protected SiteModel.Base m_siteModel;
    protected double[] m_branchLengths;
    protected double originTime;
    protected int nodeCount;


    public Hashtable<Integer,List<List<Integer>>> ancestralStates ;
    public double[][] partialLikelihoods ;
    protected double[] scalingFactors;
    protected boolean useScaling = false;


    private double scalingThreshold = 1.0E-100;



    @Override
    public void initAndValidate() {

        nodeCount = treeInput.get().getNodeCount();
        m_siteModel = (SiteModel.Base) siteModelInput.get();
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
        if (originTimeInput.get() != null) {
            originTime = originTimeInput.get().getValue();
        }
        //TODO clean up
        else {
            originTime = 0.0;
        }

        //TODO rename to partial likelihoods
        partialLikelihoods = new double[nodeCount][];
        scalingFactors = new double[nodeCount];



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
        traverseLikelihood(tree.getRoot());



        if(originTime == 0.0) {
            //sum of all partial likelihoods at the root
            logP = Math.log(Arrays.stream(partialLikelihoods[tree.getRoot().getNr()]).sum()) + getLogScalingFactor();
            return logP;
        }
        else {
            //the tree log likelihood is the log(p) of unedited state at the origin
            //TODO check that origin logP is calculated
            logP = Math.log(calculateOriginPartial(tree.getRoot())) + getLogScalingFactor();
            return logP;

        }

    }

    protected void scalePartials(int nodeNumber) {

            double scaleFactor = 0.0;

            //find the highest partial likelihodo
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


    public void traverseAncestral(Node node) {

            if (node.isLeaf() ) {
                List<List<Integer>> possibleLeafAncestors = get_possible_ancestors(dataInput.get().getCounts().get(node.getNr()));
                ancestralStates.put(node.getNr(), possibleLeafAncestors);

            } else {

                final Node child1 = node.getLeft();
                final Node child2 = node.getRight();

                traverseAncestral(child1);
                traverseAncestral(child2);

                List<List<Integer>> AncSetChild1 = ancestralStates.get(child1.getNr());
                List<List<Integer>> AncSetChild2 = ancestralStates.get(child2.getNr());


                List<List<Integer>> AncSetNode = new ArrayList<>(AncSetChild1);
                // intersection of children ancestral states
                AncSetNode.retainAll(AncSetChild2);

                ancestralStates.put(node.getNr(), AncSetNode);

            }


    }

    public void traverseLikelihood(Node node) {

        //TODO remove this because not needed (origin handled separately now)
        if( node != null ) {

            if (node.isLeaf()) {
                //TODO partials
                double[] LeafpartialLikelihoods = init_partialLikelihoods_leaf(ancestralStates.get(node.getNr()).size());
                partialLikelihoods[node.getNr()] = LeafpartialLikelihoods;


            } else {

                final Node child1 = node.getLeft();
                final Node child2 = node.getRight();

                traverseLikelihood(child1);
                traverseLikelihood(child2);

                double[] partials = calculatePartials(node.getNr(),child1,child2);
                partialLikelihoods[node.getNr()] = partials;

                if (useScaling) {
                    scalePartials(node.getNr());
                }

            }
        }



    }

    public double[] calculatePartials(int nodeNr, Node child1, Node child2) {
        //here. implement felsensteins's pruning algorithm

        //initialize an array for the partials
        double[] partials = new double[ancestralStates.get(nodeNr).size()];

        //todo potentially remove the following dead code: when a List<List<Integer>> contains a single List<Integer> it is not considered a List<List>>

//        if (ancestralStates.get(nodeNr).size() ==5) {
//            //root node
//             List<Integer> start_state = new ArrayList<Integer>() {{
//                add(0);
//                add(0);
//                add(0);
//                add(0);
//                add(0);
//            }}; ;
//            double child1partialsum = sum_partial_child(start_state, child1Nr);
//            double child2partialsum = sum_partial_child(start_state, child2Nr);
//
//            double product = child1partialsum * child2partialsum;
//            partials[0] = product;
//
//        }
//        else {
            for (int state_index = 0; state_index < ancestralStates.get(nodeNr).size(); ++state_index) {
                List<Integer> start_state = ancestralStates.get(nodeNr).get(state_index);
                double child1partialsum = sum_partial_child(start_state, child1);
                double child2partialsum = sum_partial_child(start_state, child2);

                double statePartialLikelihood = child1partialsum * child2partialsum;
                partials[state_index] = statePartialLikelihood;
            }
//        }
        return partials;

    }

    public double calculateOriginPartial(Node rootNode) {

        //on the root node!
        List<Integer> start_state = Arrays.asList(0,0,0,0,0);
        double partialAtOrigin = sum_partial_child(start_state, rootNode);

        return partialAtOrigin;

    }

    //TODO get partial likelihood per state
    public double sum_partial_child(List<Integer> start_state, Node childNode) {
        final double branchRate = branchRateModel.getRateForBranch(childNode);

        double distance;

        //init evolutionary distance
        if (childNode.isRoot()) {
            double stemLength = originTime - childNode.getHeight();
            // TODO rename distance
            distance = stemLength * branchRate ;

        }
        else {
            distance = childNode.getLength() * branchRate ;
        }

        // calculate partials
        // TODO rename likelihood to transition from the start state into any end state at the child node?
        // Name: StatePartialLikelihood
        double sum = 0;
        if(childNode.isLeaf()) {

            List<Integer> end_state = ancestralStates.get(childNode.getNr()).get(0);
            //TODO erase sum +
            sum = sum + substitutionModel.getSequenceTransitionProbability(start_state, end_state, distance);

        }
        else {

            for (int end_state_index = 0; end_state_index < ancestralStates.get(childNode.getNr()).size(); ++end_state_index) {

                List<Integer> end_state = ancestralStates.get(childNode.getNr()).get(end_state_index);

                // if the end state has non null partial likelihood
                if (! (partialLikelihoods[childNode.getNr()][end_state_index] == 0.0)) {

                    sum = sum + substitutionModel.getSequenceTransitionProbability(start_state, end_state, distance) *
                            partialLikelihoods[childNode.getNr()][end_state_index];
                }
            }
        }

        return sum;
    }

    //TODO CamelCase
    //TODO rename leaf partial likelihoods
   public double[] init_partialLikelihoods_leaf(int size) {

        double[] leafPartials = new double[size];
        leafPartials[0] = 1;
        return leafPartials;
   }


   //TODO CamelCase
    public static List<List<Integer>> get_possible_ancestors(List<Integer> sequence) {
        // to get all possible ancestors we just remove edits 1 by 1 along the barcode, starting from the last edited position
        List<List<Integer>> ancestors = new ArrayList();
        //adding the sequence itself as an ancestor
        ancestors.add(sequence);
        List<Integer> ancestor = new ArrayList<>(sequence);
        for(int i = sequence.size()-1; i >= 0; --i) {
            if(sequence.get(i) != 0) {
                //append possible ancestor list
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
